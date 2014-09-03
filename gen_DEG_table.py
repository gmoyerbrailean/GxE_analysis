import sys

## TODO: use sys.argv to get base folder
degBase="/wsu/home/groups/piquelab/charvey/GxE/differential_expression/DEseq2_results/"

## TODO: use sys.argv to get plate list
plateList = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12]

outTable = ''
## Get DEG data from shallow sequencing
for p in plateList:
	try:
		fd = open(degBase + "out_data_P%s/data_DEG/P%s_DEG.txt" % (p, p))
		dat = fd.read().split('\n')
		if '' in dat:
			dat.remove('')
		for l in dat[1:]: # skip header
			l = l + '\tP%s\tshallow\n' % p
			outTable += l
		fd.close()
	except:
		continue
for p in plateList:
	try:
		fd = open(degBase + "out_data_DP%s/data_DEG/DP%s_DEG.txt" % (p, p))
		dat = fd.read().split('\n')
		if '' in dat:
			dat.remove('')
		for l in dat[1:]: # skip header
			l = l + '\tP%s\tdeep\n' % p
			outTable += l
		fd.close()
	except:
		continue

## output DEG table
fd = open('DEG_table.txt', 'w')
fd.write('Treatment.ID\tTreatment\tFDR.01\tFDR.05\tFDR.1\tFDR.2\tpi0\tPlate\tDepth\n')
fd.write(outTable)
fd.close()