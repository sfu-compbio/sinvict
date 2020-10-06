import sys

if len(sys.argv) != 3:
	print("Usage:\t" + sys.argv[0] + "\t<input_sinvict_file_path>\t<output_filename>")
	exit(0)

variants = {}
with open(sys.argv[1]) as infile:
	for line in infile:
		line = line.rstrip()
		tokens = line.split()

		chromosome = tokens[0]
		position = tokens[1]
		ref = tokens[3]
		alt = tokens[5]

		key = chromosome + ":" + position
		if key not in variants:
			variants[key] = ref + "\t" + alt

outfile = open(sys.argv[2], "w")
outfile.write("##fileformat=VCFv4.3\n")
outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

for key, value in variants.items():
	key_tokens = key.split(':')
	outfile.write(key_tokens[0] + "\t" + key_tokens[1] + "\t" + "." + "\t" + value + "\t" + "." + "\t" + "PASS" + "\t" + "." + "\n")
