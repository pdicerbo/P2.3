#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void bit_set(int* a, int i){

  *a |= 1 << i;
  
}

int bit_check(int* a, int i){

  return *a >> i & 1;

}

void bit_shift(int* a, int i){
  (*a) <<= i;
}

void bit_clean(int* a, int i){

  (*a) &= ~(1 << i);
  
}

int extract_bits(int a, int start, int count){

  return (a >> start) | ((1 << count) - 1);
}

int main(){
  int j, test = 0, new_test;

  bit_set(&test, 3);

  printf("\n\t%d\n", bit_check(&test, 3));
  test = 0;
  for(j = 0; j < 3; j++)
    bit_set(&test, j);

  printf("\n\t%d\n", test);

  bit_shift(&test, 2);
  printf("\n\t%d\n", test);

  new_test = extract_bits(test, 2, 3);
  printf("\n\tnew_test = %d\n", new_test);
  
  return 0;
}
