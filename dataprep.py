a=int(input("start? "))
z=int(input("end? "))
while(a<=z):
	with open("effectfile."+str(a)) as f:
		print("read effect file "+str(a))
		pairs = f.read()
		pairs = pairs.split()		
		
	with open("msfile."+str(a)) as ms:
		print("read MS file " + str(a))
		msinfo = ms.read()
		msinfo = msinfo.split()

	def find_pos(msinfo): #finds positions and stores the snps that follow
		correct_sites = []
		snps = []
		first = msinfo.index("positions:")
		msinfo = msinfo[first+1:]
		second = msinfo.index("segsites:")	
		segsites = int(msinfo[second+1])
		start = msinfo.index("positions:")
		#this is where I expect the list of affected positions to start
		end = start + segsites
		shift=end+1
		#print(msinfo[end+1])
		while (start<=end):
			correct_sites.append(msinfo[start])
			start+=1	
		while(shift<len(msinfo)):
			snps.append(msinfo[shift])
			shift+=1	
		return (correct_sites[1:], snps, segsites)#cuts off "positions:"

	def check(segsites=find_pos(msinfo)[2], effectsfile=pairs): #checks if there is the incidence where there are more effects than the ms file lists (due to rare sites)
		return segsites==((len(effectsfile)-1)/2)

	def pair_check(effects=pairs, correct_sites=find_pos(msinfo)[0]): #to make sure that every site in the ms file is in the effects file
		try:
			for site in correct_sites:
				key=effects.index(site)
			return True
		except ValueError:
			print(str(site)+" is not in the effects file.")

	def check_ms(effects=pairs[1:], correct_sites=find_pos(msinfo)[0]): #
		try:
			i=0
			while(i<len(effects)):
				key=correct_sites.index(effects[i])
				i+=2
			return True
		except ValueError:
			print(str(effects[i])+" is an additional site in the msfile.")

	def matching(effects=pairs[1:],correct_sites=find_pos(msinfo)[0]): #this function stores the positions and their respective effects in a tuple of array for practical storage purposes. to be converted later
		positions = []
		e = []
		try:					
			i=0			
			for site in correct_sites:
				key=effects.index(site)		
				value=key+1
				positions.append(effects[key])
				e.append(effects[value])
				i+=1
			return (positions,e)
		except ValueError:
		#else:
			print("Please check effectfile"+str(a)+" for rare sites. It has "+str((len(effects))/2)+" sites instead of "+ str(find_pos(msinfo)[2]))
			realigned=realign(i)
			matching(realigned)			

	def pos_check(read=find_pos(msinfo)[0],to=pairs): #this is more of a sanity check
		store = []
		try:
			for pos in read:
				key=pairs.index(pos)
				store.append(key)
			print("From check: "+str(len(store)))
		except ValueError:
			print(str(pos)+" is not listed in the effectfile")

	def reverse_check(read=pairs[1:], to=find_pos(msinfo)[0]):
		store = []	
		try:
			i=0
			while(i<len(read)):
				key=to.index(read[i])
				store.append(key)
				i+=2
			print("From reverse check: "+str(len(store)))
			return True
		except ValueError:			
			print(str(read[i+2]))			
			print(str(read[i])+ " is the additional site in the effectfile. It is the "+ str(i/2)+'th element in the file. Realigning...')
			print(realign(i)[i])
		
	def realign(i,read=pairs[1:]):
		del read[i]
		del read[i]
		return read

	reverse_check()

	with open("effectfiletest."+str(a),"w") as f:
		try:		
			positions = matching()[0]
			#print(positions)
			effects = matching()[1]
			i = 0			
			while i < len(positions):		
				f.write(str(positions[i])+"\t \t"+str(effects[i]))
				f.write("\n")
				i+=1
		except TypeError:
			print("i stopped at"+ str(i))
			positions=matching()[0]
			#effects=matching()[1]
			effects=realign(i)
			j = 0			
			while i < len(positions):		
				f.write(str(positions[j])+"\t \t"+str(effects[j]))
				f.write("\n")
				j+=1

	with open("matrixprep"+str(a)+".txt","w") as f:	
		snps_break = find_pos(msinfo)[1]
		for snps in snps_break:
			for char in snps:
				f.write(str(char)+" ")
			f.write("\n")

	#with open("msfile."+str(a),"w") as f:

	a+=1