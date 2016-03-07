from tools import *
import numpy as np

i = 0

i = bit_set(i, 3, 1)

print(i)

t = bit_test(i, 3)

if t == True:
    print("correct")

else:
    print("wrong")
a = 0
for i in range(0,3):
    a = bit_set(a, i, 1)

print("first three bits = 1: ", a)
a = bit_shift(a,2)
print("bit_shift 2: ", a)

t = bit_test(a, 4)

if t == True:
    print("bit 4 is equal to 1")

else:
    print("bit 4 is not equal to 1")

new = extract_bits(a, 2, 3)
print(new)
print(bin(new))
new = bit_set(new, 2, 0)
print(new)
print(bin(new))

last = logical_mask(a, 16)
print(last)
print(bit_left(last, 3))


