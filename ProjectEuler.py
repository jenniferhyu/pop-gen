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
def factorial(n):
	if n<1:
		return 1
	else:
		return n*factorial(n-1)
def fib(n):
	if n==1:
		return 1
	if n==2:
		return 2	
	else:
		return fib(n-1) + fib(n-2)
def is_even(n):
	return n%2==0
def even_fib(n):
	to_map = []
	for i in range(1,n+1):
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
def largest_series(): #n is a string
	n="7316717653133062491922511967442657474235534919493496983520312774506326239578318016984801869478851843858615607891129494954595017379583319528532088055111254069874715852386305071569329096329522744304355766896648950445244523161731856403098711121722383113622298934233803081353362766142828064444866452387493035890729629049156044077239071381051585930796086670172427121883998797908792274921901699720888093776657273330010533678812202354218097512545405947522435258490771167055601360483958644670632441572215539753697817977846174064955149290862569321978468622482839722413756570560574902614079729686524145351004748216637048440319989000889524345065854122758866688116427171479924442928230863465674813919123162824586178664583591245665294765456828489128831426076900422421902267105562632111110937054421750694165896040807198403850962455444362981230987879927244284909188845801561660979191338754992005240636899125607176060588611646710940507754100225698315520005593572972571636269561882670428252483600823257530420752963450"	
	start = 0
	next = 13
	possibles = []
	#try:
	while(next<len(n) or start<(len(n)-13)):			
		substring = n[start:next]
		#print(substring)
		if "0" in substring:
			start+=1
			next=start+13
		else:
			product = 1
			for char in substring:
				product*=int(char)
			possibles.append(product)
			#print(possibles)
			start+=1
			next=start+13
	return max(possibles)
def large_sum():
	with open("large_sum.txt") as f:
		array = f.read()
		array
	return str(sum(array))[0:11]