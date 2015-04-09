RNAentropy - Clotelab - Boston College

/******************************************************************************
 *   Copyright (C) 2014  Juan Antonio Garcia Martin, Peter Clote              *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 ******************************************************************************/

Program computes structural entropy (H) for a given sequence S in 2 ways:
  1. Computing expected energy <E>(S) by dynamic programing:
      Q(S) = sum_over_all_structures(boltzman_factor(s) * E(s)
      Z(S) = Partition function(S)
      <E>  = Q(S)/Z(S)
  2. Computing expected energy <E>(S) estimating d/dT ln(Z(T)) for a given interval of T
      Q(S) = Partition function(S) uncoupling formal temperature and using T+delta_T as formal temperature.
      Z(S) = Partition function(S) 
      <E>  =  RT² * ((ln(Q(S))-ln(Z(S)))/delta_T)
  Structural entropy H(s)=<E>/RT + ln(Z(S))

Usage: ./RNAentropy "sequence" -s sequence -t temperature -d delta_Temp (compute structural entropy using <E> = RT² * d/dT ln(Z(T)) ) -z energy_is_zero [0|1] (default 0) -v (verbose, extended output)
       ./RNAentropy -h (detailed help)  

Input paramters are:
   <sequence>         : If FIRST argument is a valid nucleotide sequence it will be used input
   -s <sequence>      : Alternative flag for input sequence
   -t <temperature>   : Temperature in ºC
   -d <delta_T>       : Temperature variation used for estimating d/dT ln(Z(T))
                        NOTE: If this parameter is provided, H is computed estimating <E> = RT² * d/dT ln(Z(T))
   -v                 : Output includes method for computing H and the name of the ouput parameters
   -z <energy_is_zero>: [0|1] If value is 1, energies are set to 0, ouput is structural entropy for the uniform case 
   -h                 : Print this message 

Output:
  - H     : Structural entropy
  - <E>   : Expected energy
  - Length: Sequence length
  - G     : Ensemble free energy (-RT*ln(Z))  

Output format is:
   H: Structural entropy 	<E>: Expected energy	Length: Sequence length	G: Ensemble free energy (-RT*ln(Z))

