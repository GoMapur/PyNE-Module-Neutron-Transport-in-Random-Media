import csv
import itertools

var = "extra"

file1 = open(var + ".csv", "rb")
file2 = open("/Users/clayshieh/Documents/MATLAB/" + var + "_m.txt", "rb")

r1 = csv.reader(file1)
r2 = csv.reader(file2)
row_count1 = sum(1 for row in r1)
row_count2 = sum(1 for row in r2)
if row_count1 != row_count2:
	print "Actual has " + str(row_count1) + " rows Expected has " + str(row_count2) + "rows."
file1.close()
file2.close()
file1 = open(var + ".csv", "rb")
file2 = open("/Users/clayshieh/Documents/MATLAB/" + var + "_m.txt", "rb")
reader1 = csv.reader(file1)
reader2 = csv.reader(file2)


row = 0
for lhs, rhs in zip(reader1, reader2):
	column = 0
	for lh, rh in zip(lhs, rhs):
		if abs(float(lh) - float(rh)) >= 1e-10:
			print "Actual: " + str(lh) + " Expected: " + str(rh) + "   at: (" + str(row) + "," + str(column) + ")"
		column += 1

	row += 1
print row