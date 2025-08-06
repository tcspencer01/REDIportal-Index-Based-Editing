#!/usr/bin/python
import sys
import pysam  ; print("Package Version pysam:",  pysam.__version__)
from collections import defaultdict
import argparse


def read_redi_file(redi_file):

	REDIportal_sites_dict = defaultdict(dict)

	with open(redi_file) as f:
		for line in f.readlines():
			bed_chrom, pos0, ref, alt = line.split()
			#print("Adding REDIportal site at", bed_chrom, pos0, ref, alt)
			REDIportal_sites_dict[bed_chrom][int(pos0)-1] = ref + ".to." + alt #Convert to 0-based indexing

	return REDIportal_sites_dict


def parse_regions_in_bam(regions_file, bam_file, REDIportal_sites, coverage_threshold, chrom):

	editing_indexes = {}

	with open(regions_file) as f:
		lines = f.readlines()
		for line in lines:
			bed_chrom, pos0, posn, name = line.split()

			if bed_chrom != chrom:
				continue

			strand = name.split(":")[2]

			print("\nLoading new dsRNA region")
			print("chr:", bed_chrom, "start_pos:", pos0, "posn:", posn, "name:", name, "strand:", strand)

			region_info = calculate_index_based_editing_for_region(bam_file, bed_chrom, int(pos0)-1, int(posn)-1, REDIportal_sites, strand, coverage_threshold) #Convert to 0-based indexing
			editing_indexes[name] = region_info
			print("Finished with region", name)

	return editing_indexes


#Helpful link: https://pysam.readthedocs.io/en/latest/api.html
def calculate_index_based_editing_for_region(bam_file, chrom, region_start, region_end, REDIportal_sites, strand, coverage_threshold):

	import pysam

	samfile = pysam.AlignmentFile(bam_file, "rb" )

	reads_counter = 0
	REDIportal_sites_counter = 0
	covered_REDIportal_sites = 0
	REDIportal_sites_coverage = 0
	total_positions_in_region = region_end - region_start

	mismatch_counter = 0
	match_counter = 0

	AA_counter = 0
	AG_counter = 0
	TT_counter = 0
	TC_counter = 0

	print("Region start:", region_start, "--> Region end:", region_end)

	for pileupcolumn in samfile.pileup(chrom, region_start, region_end, truncate=True):

		reads_counter += pileupcolumn.n
		#print("\ncoverage at region %s = %s" % (pileupcolumn.pos, pileupcolumn.n))

		if int(pileupcolumn.pos) not in REDIportal_sites[chrom]: #Double-check 0 vs 1-based indexing here
			continue

		REDIportal_sites_coverage_temp = 0

		#LEAVE THE LINE BELOW COMMENTED OUT - we'll just filter later
		#if pileupcolumn.n < 6:
			#continue
		
		REDIportal_sites_counter += 1


		for pileupread in pileupcolumn.pileups: #Going through each read
			
			if not pileupread.is_del:

				REDIportal_sites_coverage_temp += 1

				mismatch_type = REDIportal_sites[chrom][int(pileupcolumn.pos)]

				#print("pileupcolumn.pos:", pileupcolumn.pos)
				#print("Found a REDIportal site with coverage", str(pileupcolumn.n))
				print("Found non-del REDIportal site", chrom+":"+str(pileupcolumn.pos)+":"+mismatch_type)

				read_base = (pileupread.alignment.query_sequence[pileupread.query_position]).upper()
				print("base:", read_base)

				if mismatch_type == "A.to.G":
					if read_base == "A":
						AA_counter += 1
					elif read_base == "G":
						AG_counter += 1
				elif mismatch_type == "T.to.C":
					if read_base == "T":
						TT_counter += 1
					elif read_base == "C":
						TC_counter += 1
				else: #In which case the editing type isn't given
					if read_base == "A":
						AA_counter += 1
					elif read_base == "G":
						AG_counter += 1
					elif read_base == "T":
						TT_counter += 1
					elif read_base == "C":
						TC_counter += 1
				

		if REDIportal_sites_coverage_temp > 0:
			covered_REDIportal_sites += 1
			REDIportal_sites_coverage += REDIportal_sites_coverage_temp

	samfile.close()

	output_strand = "."

	if strand == "+":
		print("Picking plus strand based on region annotation")
		match_counter += AA_counter
		mismatch_counter += AG_counter
		output_strand = "+"
	elif strand == "-":
		print("Picking minus strand based on region annotation")
		match_counter += TT_counter
		mismatch_counter += TC_counter
		output_strand = "-"
	elif (AA_counter + AG_counter) > (TT_counter + TC_counter):
		print("Picking plus strand based on read counts")
		match_counter += AA_counter
		mismatch_counter += AG_counter
		output_strand = "+"
	elif (TT_counter + TC_counter) > (AA_counter + AG_counter):
		print("Picking minus strand based on read counts")
		match_counter += TT_counter
		mismatch_counter += TC_counter
		output_strand = "-"

	print("AA_counter:", AA_counter)
	print("AG_counter:", AG_counter)
	print("TT_counter:", TT_counter)
	print("TC_counter:", TC_counter)

	if REDIportal_sites_counter == 1:
		print("Only one REDI site in region", chrom+":"+str(region_start+1)+":"+str(region_end+1))
	if REDIportal_sites_counter > 1:
		print("Found more than one REDI site in region", chrom+":"+str(region_start+1)+":"+str(region_end+1))

	if REDIportal_sites_coverage < coverage_threshold:
		print("Region", chrom+":"+str(region_start+1)+":"+str(region_end+1), "eliminated due to low coverage of", str(REDIportal_sites_coverage))
		
		if REDIportal_sites_counter == 0:
			print("Note that there were no REDI sites found covered in this region")
		if mismatch_counter + match_counter == 0:
			print("Note that in this region, mismatch_counter + match_counter == 0")

		return None

	if REDIportal_sites_counter == 0:
		print("No REDI sites in region", chrom+":"+str(region_start+1)+":"+str(region_end+1))
		return None

	if mismatch_counter + match_counter == 0:
		print("mismatch_counter + match_counter == 0")
		return None

	print("match_counter =", str(match_counter))
	print("mismatch_counter =", str(mismatch_counter))
	print("Final count:", str(mismatch_counter/(mismatch_counter + match_counter)))

	return [float((mismatch_counter/(mismatch_counter + match_counter))), int(REDIportal_sites_counter), int(covered_REDIportal_sites), int(REDIportal_sites_coverage), int(reads_counter), int(total_positions_in_region), str(output_strand)]

def write_output_file(editing_indexes, bam_file_path, output_bam_suffix, chrom):

	fh = open(bam_file_path+"."+chrom+output_bam_suffix, 'w+')

	fh.write("region_name"+"\t"+"editing_index"+"\t"+"total_REDI_sites"+"\t"+"covered_REDI_sites"+"\t"+"REDI_sites_total_coverage"+"\t"+"total_reads"+"\t"+"total_positions_in_region"+"\t"+"output_strand"+"\n")

	for region in editing_indexes:
		if editing_indexes[region] == None:
			continue

		fh.write(region)

		for info in editing_indexes[region]:
			fh.write("\t")
			fh.write(str(info))

		fh.write("\n")
		
	fh.close()

	print("job completed")


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description = 'calculates editing at REDIportal sites')

	parser.add_argument("-c", "--chr", dest = "input_chr",
			help = "chromosome to focus on", default = None)

	parser.add_argument("-b", "--input_bam", dest = "input_bam",
			help = "input bam file", default = None)

	parser.add_argument("-o", "--output_suffix", dest = "input_suffix",
			help = "output bam file Suffix", default = ".editing_levels.txt")

	parser.add_argument("-r", "--rediportal", dest = "redi_sites",
			help = "list of REDIportal sites in BED format", default = None)

	parser.add_argument("-g", "--genomic_regions", dest = "genomic_regions",
		help = "file containing genomic regions for which the script will calculate index-based RNA editing levels", default = None)

	parser.add_argument("-t", "--coverage_threshold", dest="input_coverage_threshold",
			help = "coverage threshold (corresponds to REDIportal_sites_coverage in this script)", default = 0) #Note that the default is 0 so the script still works if you forget to set this

	input_args = parser.parse_args()

	chrom=str(input_args.input_chr)
	bam_file_path=str(input_args.input_bam)
	REDI_sites_path=str(input_args.redi_sites)
	regions_file=str(input_args.genomic_regions)
	coverage_threshold=int(input_args.input_coverage_threshold)
	output_bam_suffix=str(input_args.input_suffix)

	REDIportal_sites = read_redi_file(REDI_sites_path)
	editing_indexes = parse_regions_in_bam(regions_file, bam_file_path, REDIportal_sites, coverage_threshold, chrom)
	write_output_file(editing_indexes, bam_file_path, output_bam_suffix, chrom)
