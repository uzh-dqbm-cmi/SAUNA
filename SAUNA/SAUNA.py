#!/usr/bin/env python
# coding: utf-8

# This code contains code from from the NucPosSimulator project by R. Schöpflin, V. B. Teif, O. Müller, C. Weinberg, K. Rippe and G. Wedemann , licensed under GPL-3.0 (https://opensource.org/license/gpl-3-0)
# NucPosSimulator is Copyright (C) 2013  Robert Schoepflin, Gero Wedemann     Contact: robert.schoepflin@fh-stralsund.de or gero.wedemann@fh-stralsund.de
# For the paper of the project see Robert Schöpflin, Vladimir B. Teif, Oliver Müller, Christin Weinberg, Karsten Rippe, Gero Wedemann, Modeling nucleosome position distributions from experimental nucleosome positioning maps, Bioinformatics, Volume 29, Issue 19, October 2013, Pages 2380–2386, https://doi.org/10.1093/bioinformatics/btt404
#The code was re-written in python and refactored 
#- to optimize memory usage and speed: 
#  1. initialization, 
#  2. use of memory mapping, 
#  3. vectorized operations instead of for loops
#  These changes enable the analysis of whole genomes, giving utility in whole-genome cfDNA-sequencing
#- to change certain parameters to accomodate the optimization steps
#- to fit python logic (instead of C++)
#- for legibility
#- new functions were added to extend the tool's utility


import types
import gc
import tempfile
import os.path as path
import sys
import os
import pandas as pd
import ctypes
import shutil
import inspect
import random
import math
from typing import List, Iterator, Optional, Tuple, Union
import numpy as np
import gzip

EPS = 2.220446049250313e-16  # Equivalent to numeric_limits<double>::epsilon()
K_B = 8.314513e-3  # in kJ/(mol * K)   GROMACS units
REFERENCE_TEMPERATURE = 293.0  # K
GENERIC_NUC_LENGTH = 147  # bp
MIN_PROBABILITY = 0.000000001

# Move constants
MAX_NUC_SHIFT = 60  # bp
MAX_NUC_PAIR_SHIFT = 60  # bp

ADD_RATE = 4*10**(-6) 
DELETE_RATE = 4*10**(-6)
SHIFT_RATE = 0.55 - 4*10**(-6) 
PAIR_SHIFT_RATE = 0.45-4*10**(-6)

# IO constants
MAX_LOCUS_LENGTH = 10_000_000_000  # bp


def process_file(file_location):
    # Check if the file exists
    if not os.path.exists(file_location):
        print("File not found!")
        return
    
    df = pd.read_csv(file_location, sep='\t', usecols=[1], header = None)
    df = df.values.flatten()
    return df


class AbstractException(BaseException):
    def __init__(self, msg, file, line):
        self.msg = msg
        self.file = file
        self.line = line

    def getMessage(self):
        return self.msg

    def getFile(self):
        return self.file

    def getLine(self):
        return self.line

class NucPosRunTimeException(AbstractException):
    def __init__(self, msg, file, line):
        super().__init__(msg, file, line)

    def __del__(self):
        pass  # No special cleanup needed in Python


class NucPosIOException(AbstractException):
    def __init__(self, msg, file, line):
        super().__init__(msg, file, line)

    def __str__(self):
        return f"NucPosIOException: {self.msg} at {self.file}:{self.line}"

    def __repr__(self):
        return f"NucPosIOException('{self.msg}', '{self.file}', {self.line})"


N = 624
M = 397
MATRIX_A = 0x9908b0df   # constant vector a
UPPER_MASK = 0x80000000  # most significant w-r bits
LOWER_MASK = 0x7fffffff  # least significant r bits

mt = (ctypes.c_uint32 * N)()
mti = N + 1  # mti==N+1 means mt[N] is not initialized

def init_genrand(s):
    global mt, mti
    mt[0] = s & 0xffffffff
    for i in range(1, N):
        mt[i] = (1812433253 * (mt[i-1] ^ (mt[i-1] >> 30)) + i) & 0xffffffff
    mti = N

def init_by_array(init_key, key_length):
    global mt, mti
    init_genrand(19650218)
    i, j = 1, 0
    k = max(N, key_length)
    while k:
        mt[i] = ((mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525)) +
                 init_key[j] + j) & 0xffffffff
        i += 1
        j += 1
        if i >= N:
            mt[0] = mt[N-1]
            i = 1
        if j >= key_length:
            j = 0
        k -= 1
    for k in range(N-1, 0, -1):
        mt[i] = ((mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941)) - i) & 0xffffffff
        i += 1
        if i >= N:
            mt[0] = mt[N-1]
            i = 1
    mt[0] = 0x80000000  # MSB is 1; assuring non-zero initial array

def genrand_int32():
    global mt, mti
    mag01 = [0x0, MATRIX_A]
    if mti >= N:
        if mti == N + 1:
            init_genrand(5489)
        for kk in range(N-M):
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK)
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1]
        for kk in range(N-M, N-1):
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK)
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1]
        y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK)
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1]
        mti = 0
    y = mt[mti]
    mti += 1
    y ^= (y >> 11)
    y ^= (y << 7) & 0x9d2c5680
    y ^= (y << 15) & 0xefc60000
    y ^= (y >> 18)
    return y

def genrand_int31():
    return genrand_int32() >> 1

def genrand_real1():
    return genrand_int32() * (1.0 / 4294967295.0)

def genrand_real2():
    return genrand_int32() * (1.0 / 4294967296.0)

def genrand_real3():
    return ((genrand_int32() >> 1) + 0.5) * (1.0 / 4294967296.0)

def genrand_res53():
    a = genrand_int32() >> 5
    b = genrand_int32() >> 6
    return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0)


class Interval:
     # Define the TYPE enumeration
    class TYPE:
        DNA = "DNA"
        NUC = "NUC"
        
    def __init__(self, begin: int, end: int, type_: str):
        assert end >= begin, "End must be greater than or equal to begin"
        self.begin = begin
        self.end = end
        self.length = end - begin
        self.type_ = type_

    def getBegin(self) -> int:
        return self.begin

    def getEnd(self) -> int:
        return self.end

    def getLength(self) -> int:
        return self.length

    def getType(self) -> str:
        return self.type_

    def isInInterval(self, position: int) -> bool:
        return self.begin <= position < self.end

    def setValues(self, begin: int, end: int):
        assert end > begin, "End must be greater than begin"
        self.begin = begin
        self.end = end
        self.length = end - begin

class Configuration:
    def __init__(self, filename: str, start_nucs_input: np.ndarray, minVal: int, length: int,locusBegin: int, locusLength: int, nucLength: int, chromosome: str):
        self.nucLength = nucLength
        self.chromosome = chromosome
        
        assert nucLength >= GENERIC_NUC_LENGTH  # bp
        assert locusBegin >= 0
        assert locusLength > 0

        self.locusBegin = locusBegin
        self.locusLength = locusLength
        
        start_nucs_input = start_nucs_input-minVal
        start_nucs_input += 147
        start_nucs_input = filter_peak_positions(start_nucs_input, 147)
        self.intervals,self.numOfNucleosomes = get_start_Intervals_new(nucLength,locusLength, start_nucs_input)

        del start_nucs_input
        gc.collect()

        self.energy = 0
        self.step = 0
        self.temperature = 0
        
    def getInterval(self, position):
        index = self.getIntervalIter(position)
        return self.intervals[index]
    
    def getIntervalIter(self, position):
        low = 0
        high = len(self.intervals) - 1
        result = high  # Initialize result to the last interval index

        while low <= high:
            mid = (low + high) // 2
            interval = self.intervals[mid]

            if interval.isInInterval(position):
                result = mid
                break
            elif interval.getEnd() <= position:
                low = mid + 1
            else:
                high = mid - 1

        assert self.intervals[result].isInInterval(position)  # Ensure that an interval was found
        return result
   
    def addNucleosome(self, startPosition: int):
        it = self.getIntervalIter(startPosition)

        assert self.intervals[it].getType() == Interval.TYPE.DNA
        assert self.intervals[it].getEnd() > startPosition + self.nucLength
        assert self.intervals[it].isInInterval(startPosition)

        begin = self.intervals[it].getBegin()
        end = self.intervals[it].getEnd()

        b1 = begin
        e1 = startPosition
        b2 = startPosition
        e2 = startPosition + self.nucLength
        b3 = startPosition + self.nucLength
        e3 = end

        self.intervals[it].setValues(b1, e1)  # Set new boundaries of DNA interval
        it += 1  # Step forward
        self.intervals = np.insert(self.intervals, it, Interval(b2, e2, Interval.TYPE.NUC), axis=0)

        self.numOfNucleosomes += 1
        self.intervals = np.insert(self.intervals, it+1, Interval(b3, e3, Interval.TYPE.DNA), axis=0)        
        
        assert b1 < e1
        assert e1 == b2
        assert b2 < e2
        assert e2 == b3
        assert b3 < e3
        
     
    def deleteNucleosome(self, index: int):
        assert 0 <= index < self.numOfNucleosomes
        it = self.getNucleosomeIter(index)
        
        assert it != 0  # Ensure it's not the first element
        it -= 1
        assert self.intervals[it].getType() == Interval.TYPE.DNA
        dna0 = self.intervals[it]
        it += 1
        assert self.intervals[it].getType() == Interval.TYPE.NUC
        nuc = self.intervals[it]
        nuc_index = it
        it += 1
        assert self.intervals[it].getType() == Interval.TYPE.DNA
        dna1 = self.intervals[it]
        it += 1  # Move one position beyond the linker DNA

        dna0.setValues(dna0.getBegin(), dna1.getEnd())
       
        self.intervals = np.delete(self.intervals, slice(nuc_index, it), axis=0)
        
        self.numOfNucleosomes -= 1
    
    def getNucleosomeIter(self, nucIndex: int):
        assert 0 <= nucIndex < self.numOfNucleosomes
        index = int((nucIndex+1)*2 -1)

        # Ensure that the nth nucleosome is found
        assert index != -1, f"Nucleosome with index {nucIndex} not found"
        return index
    
    def shiftNucleosome(self, index: int, distance: int):
        assert 0 <= index < self.numOfNucleosomes
        if distance == 0:
            return

        it = self.getNucleosomeIter(index)
        
        assert it != 0  # Ensure it's not the first element
        it -= 1
        assert self.intervals[it].getType() == Interval.TYPE.DNA
        dna0 = self.intervals[it]

        it += 1
        assert self.intervals[it].getType() == Interval.TYPE.NUC
        nuc = self.intervals[it]

        assert it != len(self.intervals)-1  # Ensure it's not the end of the list
        it += 1
        assert self.intervals[it].getType() == Interval.TYPE.DNA
        dna1 = self.intervals[it]

        # Shift to the left
        if distance < 0:
            assert dna0.getEnd() - dna0.getBegin() > (-1 * distance)
           
        else:  # Shift to the right
            assert dna1.getEnd() - dna1.getBegin() > distance
           

        b1 = dna0.getBegin()
        e1 = dna0.getEnd() + distance
        b2 = nuc.getBegin() + distance
        e2 = nuc.getEnd() + distance
        b3 = dna1.getBegin() + distance
        e3 = dna1.getEnd()
        dna0.setValues(b1, e1)
        nuc.setValues(b2, e2)
        dna1.setValues(b3, e3)
        

    def getNucleosomeInterval(self, nucIndex: int) -> Interval:
        it = self.getNucleosomeIter(nucIndex)
        return self.intervals[it]

    def getStartPositionOfNuc(self, nucIndex: int) -> int:
        assert 0 <= nucIndex < self.numOfNucleosomes
        index = self.getNucleosomeIter(nucIndex)
        it = self.intervals[index]
        return it.getBegin()
        
    def isStartPositionFree(self, pos: int) -> bool:
        end = self.intervals[-1].getEnd()
      
        result = False

        assert 0 <= pos < end
        it = self.getIntervalIter(pos)
        interval = self.getInterval(it) 

        if interval.getType() == Interval.TYPE.DNA and interval.getBegin() < pos and interval.getEnd() - self.nucLength > pos:
            result = True
        return result
       

    def canShiftNucleosome(self, nucIndex: int, distance: int) -> bool:
        assert 0 <= nucIndex < self.numOfNucleosomes

        result = True
        it = self.getNucleosomeIter(nucIndex)
        assert self.intervals[it].getType() == Interval.TYPE.NUC
        nuc = self.intervals[it]

        if distance < 0:  # Test left shift
            if it == 0:
                result = False
            else: 
                it -= 1
                assert self.intervals[it].getType() == Interval.TYPE.DNA
                dna = self.intervals[it]
                if nuc.getBegin() + distance <= dna.getBegin():
                    result = False
        elif distance > 0:  # Test right shift
            if it == len(self.intervals)-1:
                result=False
            else:
                it += 1
             
                assert self.intervals[it].getType() == Interval.TYPE.DNA
                dna = self.intervals[it]
                if nuc.getEnd() + distance >= dna.getEnd():
                    result = False

        return result
        
    def getNumOfFreePositions(self) -> int:
        freePositions = 0
        for interval in self.intervals:
            if interval.getType() == Interval.TYPE.DNA and interval.getLength() >= self.nucLength + 2:
                freePositions += interval.getLength() - self.nucLength - 2
        return freePositions

    def canShiftNucleosomePair(self, nucIndex0: int, nucIndex1: int, distance: int) -> bool:
        assert (nucIndex0 == nucIndex1 - 1) or (nucIndex0 == self.numOfNucleosomes - 1 and nucIndex1 == 0)
        assert 0 <= nucIndex0 < self.numOfNucleosomes
        assert 0 <= nucIndex1 < self.numOfNucleosomes
        
        result = False
        
        if nucIndex0 == nucIndex1 - 1:
            if distance < 0:
                result = self.canShiftNucleosome(nucIndex0, distance)
            elif distance >0:
                result = self.canShiftNucleosome(nucIndex1, distance)
                
        elif nucIndex0 == self.numOfNucleosomes - 1 and nucIndex1 == 0:
            result = False
       
        else:
            raise AssertionError("Mismatching indices in PairShiftMove")
        return result
    

    def getNucIndex(self, pos: int) -> int:
        assert 0 <= pos < self.intervals[-1].getEnd()
        interval_index = self.getIntervalIter(pos)
        nucindex = int((interval_index+1)/2 -1 )
        return nucindex
    
    def setStep(self, step: int):
        assert step > 0
        self.step = step
        
    def getNucLength(self) -> int:
        return self.nucLength
    
    def getNumOfNucleosomes(self) -> int:
        return self.numOfNucleosomes
    
    def setTemperature(self, temperature: float):
        self.temperature = temperature
        
    def increaseSteps(self):
        self.step += 1
    
    def addDeltaEnergy(self, deltaEnergy):
        self.energy += deltaEnergy
    
    def getLength(self) -> int:
        return self.intervals[-1].getEnd()
    
    def getEnergy(self) -> float:
        return self.energy
    
    def getStep(self):
        return self.step
    
    def getChromosome(self):
        return self.chromosome
    
    def getLocusBegin(self):
        return self.locusBegin
    
    def getTemperature(self):
        return self.temperature

def filter_peak_positions(peak_positions, min_distance):
    while True:
     
        peak_positions = peak_positions[np.concatenate(([True], np.diff(peak_positions) > min_distance))]
        if np.all(np.diff(peak_positions) > min_distance):
            break
            
    return peak_positions

def get_start_Intervals_new(nucLength,locusLength, indices):
    half_nucleosome = int(nucLength/2)
    intervals = [Interval(0,indices[0]-half_nucleosome , Interval.TYPE.DNA),Interval(indices[0]-half_nucleosome,indices[0]-half_nucleosome+ nucLength,Interval.TYPE.NUC), Interval(indices[0]-half_nucleosome+ nucLength,indices[1]-half_nucleosome,Interval.TYPE.DNA )]
    numOfNucleosomes = len(indices)
    for number,index in enumerate(indices):
        if number == 0:
            continue
        start = intervals[number*2].getEnd()
        if number != len(indices)-1:
            subintervals = [Interval(start,start + nucLength,Interval.TYPE.NUC ),Interval(start + nucLength,indices[number+1]-half_nucleosome,Interval.TYPE.DNA)]
            intervals.extend(subintervals)
        if number == len(indices)-1:
    
            if index + half_nucleosome >= locusLength:
                begin = intervals[-1][-1].getBegin()
                intervals[-1][-1].setValues(begin, locusLength)
            else:
                
                subintervals = [Interval(start, start + nucLength, Interval.TYPE.NUC), Interval(start + nucLength, locusLength,Interval.TYPE.DNA)]
                intervals.extend(subintervals)
    return np.array(intervals),numOfNucleosomes

class Energy:
    def __init__(self, parent_dir: str, filename: str, probabilities: np.ndarray, locusBegin: int, locusEnd: int, bindingEnergy: float):
        # Create NumPy array for energy values
        self.locusBegin = locusBegin
        self.locusEnd = locusEnd
        self.bindingEnergy = bindingEnergy
        def get_penalty(positions, data, distance, penalty_scale = 1):
            total_penalty = np.maximum(0, data[positions - distance] - data[positions]) + np.maximum(0, data[positions + distance] - data[positions])
            decayed_penalty = np.exp(-total_penalty * penalty_scale)
            del total_penalty
            gc.collect()
            return decayed_penalty

        def calculate_lower_proximity(positions, data, distance, penalty_scale = 1):
            total_penalty = get_penalty(positions[(positions - distance >= 0) & (positions + distance < len(data))], data, distance, penalty_scale)
            probabilities = np.zeros(len(positions))
            probabilities[(positions - distance >= 0) & (positions + distance < len(data))] = total_penalty
            del total_penalty
            gc.collect()
            return probabilities
       
        probabilities = calculate_lower_proximity(np.array(range(0,len(probabilities))), probabilities, 40,3000) * probabilities

        probabilities = np.clip(probabilities, MIN_PROBABILITY, None)
        if np.any(probabilities < MIN_PROBABILITY) or np.any(probabilities > 1.0 + EPS):
            errMsg = f"Incorrect probability value in Energy construction\n"
            errMsg += f"Value in the range [{MIN_PROBABILITY}, 1.0]"
            frame = inspect.currentframe()
            line_number = frame.f_lineno
            raise NucPosRunTimeException("I", __file__, line_number)
        
        # Calculate energy values using vectorized operations
        filename = filename +".energy.dat"
        # Create a temporary directory within the specified parent directory
        self.temp_dir = tempfile.mkdtemp(dir=parent_dir)
        filename = path.join(self.temp_dir, filename)
       
        self.energyValues = np.memmap(filename,dtype = "float32", mode='w+', shape=probabilities.shape)
        self.filename = filename
        del filename
        self.energyValues[:] = -1.0 * (np.log(probabilities) * K_B * REFERENCE_TEMPERATURE)
        assert len(probabilities) == len(self.energyValues)
        assert locusEnd - locusBegin == len(self.energyValues)
        del probabilities #remove because uses much memory
        gc.collect()
        self.energyValues.flush()
        
    def cleanup(self):
        # Remove the temporary directory and its contents
        shutil.rmtree(self.temp_dir)
        
    def getEnergy(self, index: int) -> float:
        assert 0 <= index < len(self.energyValues)
        return self.energyValues[index]

    def get_binding_energy(self) -> float:
        return self.bindingEnergy  # Simply return the binding energy attribute
    
    def getShiftEnergyDifference(self, fromCenterPos: int, toCenterPos: int) -> float:
        delEnergy = -self.getEnergy(fromCenterPos)
        addEnergy = self.getEnergy(toCenterPos)
        return delEnergy + addEnergy

    def getDeleteEnergyDifference(self, centerPos: int) -> float:
        return -self.getEnergy(centerPos) - self.bindingEnergy

    def getAddEnergyDifference(self, centerPos: int) -> float:
        return self.getEnergy(centerPos) + self.bindingEnergy
    def get_all_energy(self):
        return self.energyValues


class EnergyFactory:
    def __init__(self, parent_dir, filename,pReads, locusBegin, locusEnd):
        assert locusBegin >= 0
        assert locusBegin < locusEnd
        assert pReads is not None
        self.filename = filename
        self.parent_dir = parent_dir
        self.locusBegin = locusBegin  # Store locusBegin as an attribute
        self.locusEnd = locusEnd  # Store locusEnd as an attribute
        
        length = locusEnd - locusBegin
        self.nucCenters = np.zeros(length)
        del length
        self.nucCenters[147:len(pReads)+147] = pReads
        del pReads
        gc.collect()
    def smooth_values(self, sigma):
        # Create the Gaussian kernel
        half_kernel = 8 * np.ceil(sigma)
        kernel_indices = np.arange(-half_kernel, half_kernel + 1)
        kernel = self.gauss(kernel_indices, sigma)
        del kernel_indices
        # Pad nucCenters to handle edge cases
        self.nucCenters = np.pad(self.nucCenters, (int(half_kernel), int(half_kernel)), mode='edge')
        del half_kernel
        gc.collect()
        # Perform convolution using np.convolve
        self.nucCenters = np.convolve(self.nucCenters, kernel, mode='valid')
        return self.nucCenters

    def gauss(self, x, sigma):
        return 1.0 / (sigma * 2.0 * np.sqrt(2 * np.pi)) * np.exp(-0.5 * (x * x) / (sigma * sigma))

    def give_energy(self, sigma, binding_energy):
        # Smooth the values
        self.nucCenters = self.smooth_values(sigma)
        # Shift the data to make all values non-negative
        min_value = np.min(self.nucCenters)
        self.nucCenters = self.nucCenters + abs(min_value)
        # Determine max and sum
        maximum = np.max(self.nucCenters)
        assert maximum > 0
        # Normalize
        self.nucCenters = self.nucCenters / maximum
        assert np.all((0 <= self.nucCenters) & (self.nucCenters <= 1.0))
        return Energy(self.parent_dir,self.filename, self.nucCenters, self.locusBegin, self.locusEnd, binding_energy)

def open_file(file_path):
    columns_to_load = [3]
    arrs = []
    chunk_size = 50000
        # Read the first row to get the first element of the first and second columns
    with gzip.open(file_path, 'rt') as f:
        first_row_df = pd.read_csv(f, sep='\t', usecols=[0,1], nrows=1, header=None)
        chromosome = first_row_df.iloc[0, 0]  # Element of the first row, first column
        start = first_row_df.iloc[0, 1]  # Element of the first row, second column
        del first_row_df
        
    with gzip.open(file_path, 'rt') as f:
        for chunk in pd.read_csv(f, sep='\t', chunksize=chunk_size, usecols=columns_to_load, header=None):
            chunk[3] = chunk[3].astype(float)
            arrs.append(chunk.values)
    arrays = np.concatenate(arrs)
    del arrs
    gc.collect()
            
    with gzip.open(file_path, 'rt') as f:
        last_row_df = pd.read_csv(f, sep='\t', usecols=[1], nrows=1, header=None, skiprows= (len(arrays)-1))
        end = last_row_df.iloc[0,0]  # Last element of the second column
        del last_row_df
        
    return arrays[:,0], start ,end,chromosome

class ReadReader:
    def __init__(self, filename: str):
        self.locusBegin = 0
        self.locusEnd = 0
        min_val = float('inf')
        max_val = float('-inf')
        self.pReads,self.minValue,self.maxValue,chrom  = open_file(filename)
        
        self.chromosome = "chr" + str(chrom)
        del chrom       
        assert self.minValue >= 0
        assert self.maxValue >= 0

        # add flanking DNA
        self.locusBegin = max(0, self.minValue)
        self.locusEnd = self.maxValue + 2*GENERIC_NUC_LENGTH

        assert len(self.pReads) > 0
        assert self.locusBegin < self.locusEnd
        assert 0 <= self.locusBegin < self.locusEnd

        print(f"Imported {len(self.pReads)} nucleosome reads successfully.\n"
              f"Read range:\t{self.minValue}:{self.maxValue}\n"
              f"Locus begin:\t{self.locusBegin}\n"
              f"Locus end:\t{self.locusEnd}")  
    
    def getReads(self) -> List[Tuple[Union[str, int], int, int]]:
        return self.pReads

    def getLocusBegin(self) -> int:
        return self.locusBegin

    def getLocusEnd(self) -> int:
        return self.locusEnd

    def getLocusLength(self) -> int:
        return self.locusEnd - self.locusBegin

    def getChromosome(self) -> str:
        return self.chromosome
    def getMin(self)->int:
        return  self.minValue
    def getMax(self)->int:
        return  self.maxValue

class ConfigWriter:
    def __init__(self):
        pass

    def writeConfig2Bed(self, config, out):
        for interval in config.intervals:
            if interval.getType() == Interval.TYPE.NUC:
                out.write(f"{config.getChromosome()}\t"
                          f"{config.getLocusBegin() + interval.getBegin()-147}\t"
                          f"{config.getLocusBegin() + interval.getEnd()-147}\n")
        out.flush()
        del config
        gc.collect()


class MoveSelector:
    def __init__(self):
        self.moves = []
        self.cumulated_probabilities = []

    def __del__(self):
        for move in self.moves:
            del move

    def addMove(self, move, probability):
        self.moves.append(move)
        self.cumulated_probabilities.append(probability)
        n = len(self.cumulated_probabilities)
        # Cumulate probabilities
        if n > 1:
            self.cumulated_probabilities[-1] += self.cumulated_probabilities[-2]

    def next_(self):
        # Check that overall probability is 1
        assert math.isclose(self.cumulated_probabilities[-1], 1.0, abs_tol=EPS)

        # Get a random number between 0 and 1
        number = genrand_real1()

        # Select next move
        p_result = None
        for i, probability in enumerate(self.cumulated_probabilities):
            if number <= probability:
                p_result = self.moves[i]
                break
        assert p_result is not None
        return p_result

    def printRates(self):
        print("-----------------------------------------\nMove acceptance rates:")
        for move in self.moves:
            print(move.getName(), "\t", move.getAcceptanceRate())


class AbstractMove:
    def __init__(self, config, energy):
        self.config = config  # Configuration object
        self.energyFunction = energy  # Energy object
        self.counter = 0  # Counter for prepared moves
        self.accepted = 0  # Counter for accepted moves
        self.prepared = False  # Flag indicating if a move is prepared

    def prepareMove(self):
        pass

    def calcDeltaEnergy(self):
        # Implement calcDeltaEnergy logic here
        pass  # Placeholder for the actual implementation

    def performMove(self):
        pass

    def reset(self):
        # Reset move state
        self.prepared = False  # Set prepared flag to False

    def getName(self):
        # Implement getName logic here
        pass  # Placeholder for the actual implementation

    def getAcceptanceRate(self):
        # Calculate acceptance rate
        if self.counter == 0:
            return 0.0  # Return 0 if no moves have been prepared
        else:
            return self.accepted / self.counter  # Return acceptance rate as a float between 0 and 1

class AddMove(AbstractMove):
    def __init__(self, config, energy):
        super().__init__(config, energy)
        self.positions = config.getLength() - config.getNucLength() - 1
        assert self.positions > 0
        self.nucStartPos = 0

    def prepareMove(self):
        self.prepared = False
        self.counter += 1
        
        assert self.config.getLength() > self.config.getNucLength()
        
        randomIndex = genrand_int32() % self.positions + 1  # +1 because 0 is not allowed
        if self.config.isStartPositionFree(randomIndex):
            self.nucStartPos = randomIndex
            del randomIndex
            self.prepared = True
        return self.prepared

    def calcDeltaEnergy(self):
        assert self.prepared == True
        center = self.nucStartPos + self.config.getNucLength() // 2
        return self.energyFunction.getAddEnergyDifference(center)

    def performMove(self):
        assert self.prepared == True
        self.accepted += 1
        self.config.addNucleosome(self.nucStartPos)
        self.prepared = False

    def __del__(self):
        pass  # Destructor doesn't contain any specific cleanup

    def getName(self):
        return "AddMove"


class DeleteMove(AbstractMove):
    def __init__(self, config, energy):
        super().__init__(config, energy)
        self.nucIndex = 0
        # num of potential positions with nucleosome coverage = length - 2 (DNA margins)
        self.positions = config.getLength() - 2
        assert self.positions > 0

    def prepareMove(self):
        self.prepared = False
        self.counter += 1

        # position 0 is not allowed, has to be DNA
        randomPos = genrand_int32() %  self.positions + 1
        interval = self.config.getInterval(randomPos)
     
        if interval.getType() == Interval.TYPE.NUC:
            self.prepared = True
            self.nucIndex = self.config.getNucIndex(randomPos)
            del randomPos
        return self.prepared

    def calcDeltaEnergy(self):
        assert self.prepared == True
        center = self.config.getStartPositionOfNuc(self.nucIndex) + self.config.getNucLength() // 2
        return self.energyFunction.getAddEnergyDifference(center)


    def performMove(self):
        assert self.prepared == True
        self.accepted += 1
        self.config.deleteNucleosome(self.nucIndex)
        self.prepared = False
    def getName(self):
        return "DeleteMove"

class ShiftMove(AbstractMove):
    def __init__(self, config, energy):
        super().__init__(config, energy)
        self.nucIndex = 0
        self.distance = 0

    def reset(self):
        super().reset()
        self.distance = 0

    def prepareMove(self):
        self.prepared = False
        self.counter += 1

        num = self.config.getNumOfNucleosomes()
        if num > 0:
            # Select a random nucleosome
            self.nucIndex = genrand_int32() % num
            self.distance = (genrand_int32() % (2*MAX_NUC_SHIFT)) - MAX_NUC_SHIFT
            # Omit the zero
            if self.distance == 0:
                self.distance += 1
            assert(self.distance != 0)
            assert -MAX_NUC_SHIFT <= self.distance <= MAX_NUC_SHIFT
            if self.config.canShiftNucleosome(self.nucIndex, self.distance) == True:
                self.prepared = True

        return self.prepared

    def calcDeltaEnergy(self):
        assert self.prepared
        nuc = self.config.getNucleosomeInterval(self.nucIndex)
        from_center_pos = nuc.getBegin() + self.config.getNucLength() // 2
        to_center_pos = nuc.getBegin() + self.distance + self.config.getNucLength() // 2

        return self.energyFunction.getShiftEnergyDifference(from_center_pos, to_center_pos)
    
    
    def performMove(self):
        assert self.prepared
        self.accepted += 1
        self.config.shiftNucleosome(self.nucIndex, self.distance)
        self.prepared = False
    def getName(self):
        return "ShiftMove"
    

class PairShiftMove(AbstractMove):
    def __init__(self, config, energy):
        super().__init__(config, energy)
        self.nucIndex0 = 0
        self.nucIndex1 = 0
        self.distance = 0

    def reset(self):
        super().reset()
        self.distance = 0

    def prepareMove(self):
        self.prepared = False
        self.counter += 1

        num = self.config.getNumOfNucleosomes()
        if num >= 2:
            self.nucIndex0 = genrand_int32() % num
            self.nucIndex1 = (self.nucIndex0 + 1) % num
        
            self.distance = (genrand_int32() % (2 * MAX_NUC_PAIR_SHIFT)) - MAX_NUC_PAIR_SHIFT
            if self.distance >= 0:
                self.distance += 1
            assert -MAX_NUC_PAIR_SHIFT <= self.distance <= MAX_NUC_PAIR_SHIFT

            if self.config.canShiftNucleosomePair(self.nucIndex0, self.nucIndex1, self.distance):
                self.prepared = True

        return self.prepared

    def calcDeltaEnergy(self):
        assert self.prepared

        nuc0 = self.config.getNucleosomeInterval(self.nucIndex0)
        nuc1 = self.config.getNucleosomeInterval(self.nucIndex1)
        
        from_center_pos_0 = nuc0.getBegin() + self.config.getNucLength() // 2

        from_center_pos_1 = nuc1.getBegin() + self.config.getNucLength() // 2

        to_center_pos_0 = nuc0.getBegin() + self.distance + self.config.getNucLength() // 2

        to_center_pos_1 = nuc1.getBegin() + self.distance + self.config.getNucLength() // 2

        energy = self.energyFunction.getShiftEnergyDifference(from_center_pos_0, to_center_pos_0) + \
                 self.energyFunction.getShiftEnergyDifference(from_center_pos_1, to_center_pos_1)

        return energy


    def performMove(self):
        assert self.prepared
        self.accepted += 1

        if self.distance < 0:
            self.config.shiftNucleosome(self.nucIndex0, self.distance)
            self.config.shiftNucleosome(self.nucIndex1, self.distance)
        elif self.distance > 0:
            self.config.shiftNucleosome(self.nucIndex1, self.distance)
            self.config.shiftNucleosome(self.nucIndex0, self.distance)

        self.prepared = False
        
    def getName(self):
        return "PairShiftMove"

class SimController:
    def __init__(self, config, energyFunction):
        self.config = config

        # initialize random number generator
        seed=27
        random.seed(seed)

        print("Seed for random number generator:", seed)

        self.moveSelector = MoveSelector()
        self.moveSelector.addMove(AddMove(config, energyFunction), ADD_RATE)
        self.moveSelector.addMove(DeleteMove(config, energyFunction), DELETE_RATE)
        self.moveSelector.addMove(ShiftMove(config, energyFunction), SHIFT_RATE)
        self.moveSelector.addMove(PairShiftMove(config, energyFunction), PAIR_SHIFT_RATE)

        self.temperature = 0.0

    def run(self, steps, stepsToSave, temperature):
        self.temperature = temperature

        infoStepSize = self.computeInfoStepSize(steps)

        for i in range(steps + 1):
            self.config.setTemperature(temperature)
            if i % infoStepSize == 0:
                print("\rProgress: {:.0f} % ".format(i / steps * 100), end="")
                sys.stdout.flush()
            self.step()

        print("\rProgress: 100 % ")
        print()
        self.moveSelector.printRates()
        print("Simulation completed")

    def runAnnealing(self, steps, stepsToSave, startTemp, endTemp):
        assert startTemp > endTemp
        assert endTemp > 0.0
        assert steps > 0

        annealingSteps = 0

        if 0 < steps < 1000:
            annealingSteps = steps
        elif 1000 <= steps < 1000000:
            annealingSteps = steps // 10
        elif steps >= 1000000:
            annealingSteps = 100000
        else:
            raise Exception("Internal annealing step error")

        print("Annealing steps:", annealingSteps)

        assert 0 < annealingSteps <= 100000

        annealingFactor = pow((endTemp / startTemp), (1.0 / annealingSteps))
        annealingStepSize = steps // annealingSteps
        self.temperature = startTemp
        self.config.setTemperature(self.temperature)

        infoStepSize = self.computeInfoStepSize(steps)

        for i in range(steps + 1):

            if i % infoStepSize == 0:
                print("\rProgress: {:.0f} %\tTemperature: {:.1f} K\t #Nucs: {:6}".format(i / steps * 100, self.temperature, self.config.getNumOfNucleosomes()), end="")
                sys.stdout.flush()

            self.step()

            if i % annealingStepSize == 0 and self.temperature > endTemp:
                self.temperature *= annealingFactor
                self.config.setTemperature(self.temperature)

        print("\rProgress: 100 %\tTemperature: {:.1f} K\t #Nucs: {:6}".format(self.temperature, self.config.getNumOfNucleosomes()))
        print()
        self.moveSelector.printRates()
        print("Simulation completed")

    def step(self):
        pMove = self.moveSelector.next_()

        success = pMove.prepareMove()
        self.config.increaseSteps()

        if success:
            deltaEnergy = pMove.calcDeltaEnergy()
            p = 1.0
            if deltaEnergy <= 0:
                p = 1.0
            else:
                p = math.exp(-deltaEnergy / (K_B * self.temperature))

            randomNumber = genrand_real1()
            if randomNumber <= p:
                pMove.performMove()
                self.config.addDeltaEnergy(deltaEnergy)
        else:
            pMove.reset()

    def computeInfoStepSize(self, steps):
        infoStepSize = 1
        maxInfoSteps = 1000
        if steps / infoStepSize > maxInfoSteps:
            infoStepSize = steps // maxInfoSteps
            if steps % maxInfoSteps != 0:
                infoStepSize += 1

        return infoStepSize

def usage():
    print("\nUsage:\n\nNucPosSimulator <peak_output.tsv> <nucleosome_center_data.tsv.gz> <params.txt> [output-path]\n\n"
          "\t<peak_output.tsv>     tsv input file with peaks generated by peak calling\n"
          "\t<nucleosome_center_data.tsv.gz>     tsv.gz input file with nucleosome center data\n"
          "\t<params.txt>    parameter file\n"
          "\t[output-path]   path to an alternative output directory (optional)\n")
def main():
    try:
        if len(sys.argv) != 3 and len(sys.argv) != 4:
           usage()
           sys.exit(0)
        file_location = sys.argv[1]      
        filename = sys.argv[2]

        parameterFilename = sys.argv[3]

        outputFilebase = filename
        
        if len(sys.argv) == 5:
            outputDir = sys.argv[4]
            if os.path.isdir(outputDir):
                outputFilebase = os.path.join(outputDir, os.path.basename(filename))
            else:
               print("The output directory does not exist.")

        simSteps = None
        stepsToSave = None
        nucLength = None
        sigma = None
        annealing = None
        startTemperature = None
        endTemperature = None
        bindingEnergy = None
        temperature = None
        
        print("\n-----------------------------------------")
        readReader = ReadReader(filename)
        locusBegin = readReader.getLocusBegin()
        locusEnd = readReader.getLocusEnd()
        chrom = readReader.getChromosome()
        length = locusEnd - locusBegin
        minVal = readReader.getMin()
        maxVal = readReader.getMax()

        with open(parameterFilename, "r") as paramFile:
            for line in paramFile:
                key, value = line.strip().split("\t")
                if key == "SimSteps":
                    if value == 'default':
                        simSteps = 5*length
                    else:
                        simSteps = int(value)
                    del length
                elif key == "NucLength":
                    nucLength = int(value)
                elif key == "SmoothingSigma":
                    sigma = float(value)
                elif key == "Annealing":
                    annealing = bool(int(value))
                elif key == "StartTemperature":
                    startTemperature = float(value)
                elif key == "EndTemperature":
                    endTemperature = float(value)
                elif key == "BindingEnergy":
                    bindingEnergy = float(value)
                elif key == "Temperature":
                    temperature = float(value)

        if None in (simSteps, stepsToSave, nucLength, sigma, annealing, bindingEnergy):
            print("Missing or invalid parameter value.")

        print("\n-----------------------------------------")

        print("Simulated Annealing : ", annealing)

        name = os.path.basename(filename)
        parent_dir = os.path.dirname(filename)

        pEnergy = EnergyFactory(parent_dir,name,readReader.getReads(), locusBegin, locusEnd).give_energy(sigma, bindingEnergy)
        
        del readReader
        gc.collect()
        print("Binding energy:", pEnergy.get_binding_energy())

        if annealing:
            print("Start temperature:", startTemperature)
            print("End temperature:", endTemperature)
        
        start_nucs = process_file(file_location)
    
        config = Configuration(name,start_nucs,minVal, maxVal - minVal,locusBegin,locusEnd-locusBegin, nucLength, chrom)
        del start_nucs        

        SimController(config, pEnergy).runAnnealing(simSteps, stepsToSave, startTemperature, endTemperature)
        pEnergy.cleanup()
        del pEnergy
        gc.collect()
        
        bedFilename = outputFilebase + ".result.bed"
        with open(bedFilename, "w") as bedOut:
            config_writer = ConfigWriter()  # Create an instance of ConfigWriter
            config_writer.writeConfig2Bed(config, bedOut)  # Call the method on the instance
       

    except NucPosIOException as e:
        print("\nERROR - IO exception:")
        print(e.getMessage(), "\n")
        return -1
    except NucPosRunTimeException as e:
        print("\nERROR - Runtime exception:")
        print(e.getMessage(), "\n")
        return -1

if __name__ == "__main__":
    main()
