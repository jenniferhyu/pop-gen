#ProjectEuler
import math
import operator
from functools import reduce
def multiples(n):
	i = 0
	total = 0
	while i < n:
		if i !=0 and (i%3==0 or i%5==0):
			total+=i
			i+=1
		else:
			i+=1
	return total
def fib(n):
	if n==0:
		return 1
	if n==1:
		return 1
	if n ==2:
		return 2
	else:
		return fib(n-1) + fib(n-2)
def is_even(n):
	return n%2==0
def even_fib(n):
	to_map = []
	for i in range(n+1):
		if is_even(fib(i)):
			to_map.append(fib(i))
	return to_map

def is_prime(n):
	for i in range(2, int(math.sqrt(n)+1)):
		if n%i == 0:
			return False
	return True

def largest_prime(x):
	listed = []
	for i in range(2, int(math.sqrt(x)+1)):
		if x%i==0 and is_prime(i):
			listed.append(i)
	return max(listed)
def summed_primes(n):
	to_sum = []
	for i in range(2, int(math.sqrt(n)+1)):
		if is_prime(i):
			to_sum.append(i)
	return sum(to_sum)
def largest_product(s):
	start = 0
	s_index = 0
	listed = []
	while s_index < len(s) and s_index < 5:
		listed.append(s[s_index])
		s_index+=1
		if 0 in listed:
			break
	new=reduce(operator.mul, listed, 1)
	if new > start:
		start = new
	return start
def nth_prime(n, num=0):
	i = 0
	while i < n:
		for x in range(2, n**2):
			if is_prime(x):
				i+=1
				num = x
	return num
	
def natural_square_sum(n):
	summed = 0
	for i in range(1, n+1):
		summed += i**2
	return summed
def squared_sum(n):
	return sum(range(1,n+1))**2
def smallest_multiple(n): #n is a list
	index = len(n)//2
	while index < len(n):
		if not is_prime(n[index]):
			n.pop(index)
			index+=1
		else:
			index+=1
	return reduce(operator.mul, n, 1)
def sum_digits(n):
	if n < 10:
		return n
	else:
		all_but_last, last = split(n)
		return sum_digits(all_but_last) + last
def split(n):
	return n//10, n%10