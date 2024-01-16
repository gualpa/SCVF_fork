/**
 * @file allvars.h
 * @brief Implementación del algoritmo de quicksort.
 *
 * Copyright 2023 Andrés Nicolás Ruiz, Sebastián Rogelio Gualpa
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may
 * be used to endorse or promote products derived from this software without specific
 * prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#include "allvars.h"
#include "qsort.h"

void qsort(struct sort *a, int start, int end)
{
  struct sort temp;
  int         left,right,marker;
  double      pivot;

  if (start < end) {

     pivot = a[start].val;
     left = start;
     right = end;

     while (left < right) {
	 while (a[right].val > pivot) {
	     right--;	 
	 }	 

	 while ((left < right) && (a[left].val <= pivot)) {
	     left++;	 
	 }

	 if (left < right) {
	    temp.val = a[left].val;
	    temp.ord = a[left].ord;
            a[left].val = a[right].val;
            a[left].ord = a[right].ord;
            a[right].val = temp.val;	    
            a[right].ord = temp.ord;	    
	 }
     }

     temp.val = a[right].val;
     temp.ord = a[right].ord;
     a[right].val = a[start].val;
     a[right].ord = a[start].ord;
     a[start].val = temp.val;
     a[start].ord = temp.ord;

     marker = right;

     qsort(a, start, marker - 1);
     qsort(a, marker + 1, end);
	  
  }

}	

