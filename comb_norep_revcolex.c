/* 26.06.2013 last modification: 26.06.2013
Copyright (c) 2013 by Siegfried Koepf

This file is distributed under the terms of the GNU General Public License
version 3 as published by the Free Software Foundation.
For information on usage and redistribution and for a disclaimer of all
warranties, see the file COPYING in this distribution.

k-combinations without repetition in reverse colexicographic order

Algorithm by Siegfried Koepf

Functions:
  int gen_comb_norep_revcolex_init(unsigned char *vector, const unsigned char n, const unsigned char k)
    Test for special cases
    Initialization of vector
    Possible return values are: GEN_ERROR, GEN_EMPTY, GEN_NEXT

  int gen_comb_norep_revcolex_next(unsigned char *vector, const unsigned char k)
    Transforms current figure in vector into its successor
    Possible return values are: GEN_NEXT, GEN_TERM

Arguments:
  unsigned char *vector; //pointer to the array where the current figure is stored
  const unsigned char n; //length of alphabet
  const unsigned char k; //length of figures

Usage and restrictions:
  Arguments and elements in vector are restricted to the interval (0, 255)
  Memory allocation for vector must be provided by the calling process
  k must be <= n

Cardinality:
  n! / ((n - k)! * k!) == Binomial(n, k)
*/

#include <stdio.h>
#include <stdlib.h>
#include "_generate.h"

int gen_comb_norep_revcolex_init(unsigned char *vector, const unsigned char n, const unsigned char k)
{
int j; //index

//test for special cases
if(k > n)
 return(GEN_ERROR);

if(k == 0)
 return(GEN_EMPTY);

//initialize
vector[0] = n - k;

for(j = 1; j < k; j++)
 vector[j] = vector[j - 1] + 1;

return(GEN_NEXT);
}

int gen_comb_norep_revcolex_next(unsigned char *vector, const unsigned char k)
{
int j; //index

//easy case, decrease leftmost element
if(vector[0] > 0)
 {
 vector[0]--;
 return(GEN_NEXT);
 }

//find leftmost element to decrease
for(j = 1; j < k; j++)
 if(vector[j] > j)
  break;

//terminate if vector[k - 1] == n - k
if(j >= k)
 return(GEN_TERM);

//decrease
vector[j]--;

//set left-hand elements
while(j > 0)
 {
 vector[j - 1] = vector[j] - 1;
 j--;
 }

return(GEN_NEXT);
}

int main(void)
{
unsigned char n           = 5;    //length of alphabet
unsigned char k           = 3;    //length of figures
unsigned char *vector     = NULL; //where the current figure is stored
int           gen_result;         //return value of generation functions
unsigned int  set_counter;        //counting generated sequences
int           x;                  //iterator

//alloc memory for vector
vector = (unsigned char *)malloc(sizeof(unsigned char) * k);
if(vector == NULL)
 {
 fprintf(stderr, "error: insufficient memory\n");
 exit(EXIT_FAILURE);
 }

set_counter = 0;
printf("comb_norep_lex(%u, %u)\n", n, k);

//initialize
gen_result = gen_comb_norep_lex_init(vector, n, k);

if(gen_result == GEN_ERROR)
 {
 fprintf(stderr, "error: couldn't initialize\n");
 return(EXIT_FAILURE);
 }

if(gen_result == GEN_EMPTY)
 {
 set_counter++;
 printf("{} (%u)\n", set_counter);
 }

//generate all successors
while(gen_result == GEN_NEXT)
 {
 set_counter++;

 for(x = 0; x < k; x++)
  printf("%u ", vector[x]);

 printf("(%u)\n", set_counter);

 gen_result = gen_comb_norep_lex_next(vector, n, k);
 }

return(EXIT_SUCCESS);
}
