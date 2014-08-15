# actions copy files into samtools2
# replace "htsjdk.samtools" with "htsjdk.samtools2" in package directive
# replace "htsjdk.samtools" with "htsjdk.samtools2" in imports and references of all files on this list
# replace SAMRecord with ReadRecord
import os, sys
import shutil
#import re

def copy(old_path, new_path, verbose=False):
	if os.path.isfile(old_path):
		if verbose:
			print("Copying %s to %s" % (old_path.replace(base_path, ""), new_path.replace(base_path, "")))
		shutil.copyfile(old_path, new_path)
	else:
		if verbose:
			print("Error: %(old_path)s not found" % locals())

def replace_in_file(file_path, replace_strings, verbose=True):
	with open(file_path) as f:
		contents = f.read()

	if verbose:
		print("%(file_path)s: " % locals())
	new_contents = contents
	for old_string, new_string in replace_strings:
		new_contents = new_contents.replace(old_string, new_string)
		for i, line in enumerate(contents.split("\n")):
			if verbose and old_string in line:
				print("     Line %(i)d:\t '%(old_string)s' => '%(new_string)s'  \n \t\t\t\t\t  in line: %(line)s" % locals())
		
	with open(file_path, 'w') as f:
		f.write(new_contents)



base_path_MacOS = "/Users/weisburd/code/picard/picard-private/Picard-public/htsjdk/src/java/htsjdk/"
#base_path_linux = "/home/unix/weisburd/code/picard-private/Picard-public/htsjdk/src/java/htsjdk/"
base_path = base_path_MacOS

package1, package2 = "htsjdk.samtools", "htsjdk.samtools2"

tool_paths_relative_to_htsjdk = ["../../../../src/java/picard//analysis/CollectFastWgsMetrics.java"]
paths_relative_to_htsjdk = [
	"samtools/SAMFileReader.java",  	# needs to extend ReadRecord
    	"samtools/SamReader.java",		# needs to extend SAMReader
    	"samtools/SamReaderFactory.java",
    	"samtools/SAMRecordFactory.java",
    	"samtools/DefaultSAMRecordFactory.java",
    	"samtools/util/SamLocusIterator.java",

	#"samtools/BAMFileReader.java", 
    "samtools/SAMRecordComparator.java", 
    "samtools/SAMRecordCoordinateComparator.java",
    "samtools/SAMRecordQueryNameComparator.java",
	"samtools/SAMFileHeader.java",
	"samtools/SAMFileSource.java",
    "samtools/SAMFileTruncatedReader.java",
    "samtools/SamFileValidator.java",
    "samtools/SAMLineParser.java",
    "samtools/SAMRecordIterator.java",
    "samtools/SAMRecordSetBuilder.java",
    "samtools/SAMRecordUtil.java",
    "samtools/SAMUtils.java", # needed within SAMRecord.java
    "samtools/GenomicIndexUtil.java",

    "samtools/AlignmentBlock.java", 
    "samtools/SAMTextHeaderCodec.java",
    "samtools/TextTagCodec.java",
    "samtools/BinaryTagCodec.java",
    "samtools/TagValueAndUnsignedArrayFlag.java",
    "samtools/AbstractSAMHeaderRecord.java",
    "samtools/SamInputResource.java",
    "samtools/SAMTextReader.java",
    "samtools/SAMTextWriter.java",
 	"samtools/SAMFileWriterImpl.java",
 	"samtools/SamPairUtil.java",
 	"samtools/SAMSequenceRecord.java",
 	"samtools/SAMSequenceDictionary.java",
 	"samtools/SAMFileWriter.java",
 	"samtools/SAMSortOrderChecker.java",
 	"samtools/util/ProgressLoggerInterface.java",
 	"samtools/util/ProgressLogger.java",
 	"samtools/util/QualityEncodingDetector.java",
 	"samtools/util/SequenceUtil.java",
 	"samtools/util/WholeGenomeReferenceSequenceMask.java",
    "samtools/SAMReadGroupRecord.java",
    "samtools/SAMHeaderRecordComparator.java",
    "samtools/SAMProgramRecord.java",
    "samtools/util/SamRecordIntervalIteratorFactory.java",

	"samtools/SAMRecord.java",		# needs to extend ReadRecord
	"samtools/filter/SamRecordFilter.java",  			# needs to return ReadRecord instead of SAMRecord
	"samtools/filter/SecondaryAlignmentFilter.java",  	# used by CollectFastWGSMetrics, needs to take and return ReadRecord instead of SAMRecord
	"samtools/filter/SecondaryOrSupplementaryFilter.java",
	"samtools/filter/DuplicateReadFilter.java",
	"samtools/filter/IntervalFilter.java",
	"samtools/filter/AggregateFilter.java",
	"samtools/filter/FilteringIterator.java", 
	"samtools/util/IntervalUtil.java", 
]

dont_copy = []

# copy files into samtools2,  change package declaration
for rel_path in paths_relative_to_htsjdk:
	old_path = base_path + rel_path
	new_path = base_path + rel_path.replace("samtools/", "samtools2/")
	if not os.path.isfile(new_path): #not [x for x in dont_copy if x in rel_path]:
		copy(old_path, new_path)

	old_package_name = "htsjdk." + os.path.dirname(rel_path).replace("/", ".")
	new_package_name = old_package_name.replace("htsjdk.samtools", "htsjdk.samtools2")
	
	replace_in_file(new_path, [("package "+old_package_name+";",    "package "+new_package_name+";")])


print("------------------")
# replace strings in the new copies of the files
class_names = ["htsjdk." + rel_path.replace(".java", "").replace("/", ".") for rel_path in paths_relative_to_htsjdk]	
replacement_strings = []
for class_name in class_names:
	replacement_strings += [(class_name, class_name.replace(".samtools", ".samtools2"))]
	print((class_name, class_name.replace(".samtools", ".samtools2")))
[sys.stdout.write(x + " => " + y + "\n") for x,y in replacement_strings]

for rel_path in paths_relative_to_htsjdk + tool_paths_relative_to_htsjdk:
	#if "SAMRecord.java" not in rel_path:
	#	replacement_strings += [("SAMRecord ", "ReadRecord ")]
	replace_in_file(rel_path, replacement_strings)

	
	

	

# print occurances     
#def replace(contents):
#contents = contents.replace("package %(package)s.", "")
