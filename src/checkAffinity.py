import sys, os

if (len(sys.argv) < 3):
	print("python3 checkAffinity.py  TRAP_output, numTFs")
else:
	TRAP_file = sys.argv[1]
	numTFs = int(sys.argv[2]) + 1
	

	counter = 0
	with open(TRAP_file, 'r') as i:
		for line in i:
			line = line.strip().split('\t')
			#print(len(line))
			if len(line) != numTFs:
				counter+=1
				#print(len(line))
				print(line)
	
				break;
	print(counter)




