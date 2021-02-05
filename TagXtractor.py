#!/usr/bin/python3
#
# TagXtractor
# Version 0.0.3
# By Oddmund Nordgård (1), based on tag_to_header by Joe Hiatt, Scott Kennedy(2)
# Brendan Kohrn and Mike Schmitt(2)
# (1) Stavanger University Hospital, Stavanger, Norway
# (2) Department of Pathology, University of Washington School of Medicine, Seattle
# See licence information in separate file
#
# Isolate molecular tags, move them from within the sequenced read to the
# header region, and remove the spacer region. The ion version is adapted for
# Ion Torrent chemistry and a simplified iDES system in Stavanger.
# This version writes the resulting fastq file to stdout for piping
# It can also read from stdin if using the infile file name -
# This version (0.2) handles mistakes in the spacer
# Current version also has an option to trim the ends of the reads


#usage: tag_to_header_ion.py [-h] --infile INFILE --outprefix OUTFILE
#                            [--taglen TAGLEN] [--spacerlen SPCLEN]
#                            [--filtspacer SPACER_SEQ]
#                            [--tagstats] [--reduce]

#optional arguments:
#  -h, --help            show this help message and exit
#  --infile INFILE       Path to FASTQ input file
#  --outprefix OUTFILE   Prefix for output files. Will prepend onto file name
#  --taglen TAGLEN       Length in bases of the duplex tag sequence.[12
#  --spacerlen SPCLEN    Length in bases of the spacer sequence between duplex
#                        tag and the start of target DNA. [5]
#  --stdout              Write output to stdout instead, and report text to logfile .log

#  --filtspacer SPACER_SEQ
#                        Optional: Filter out sequences lacking the inputed
#                        spacer sequence. Not recommended due to significant
#                        base calling issues with the invariant spacer sequence
#  --tagstats            Optional: Output tagstats file and make distribution
#                        plot of tag family sizes. Requires matplotlib to be
#                        installed.
#  --reduce              Optional: Only output reads that will make a final DCS
#                        read. Will only work when the --tagstats option is
#                        invoked.


import sys
from argparse import ArgumentParser
from collections import defaultdict
from difflib import SequenceMatcher



###############################################
#FUNCTION DEFINITIONS
###############################################



# FASTQ_GENERAL_ITERATOR
# this function is a generator that returns single reads from a fastq file as
# tuplets of title, sequence and quality

def fastq_general_iterator(read1_fastq):

    #renaming of the readline function for infile
    read1_readline = read1_fastq.readline

    while True:
       #read the first line of infile
        rline = read1_readline()

        #If readline was not successful, abort function
        if not rline:
            return

        #Break while loop if it looks like a proper fastq file
        if rline[0] == '@':
                        break
        #If not a proper fastq file, check whether an int can be found and raise error
        if isinstance(rline[0], int):
            raise ValueError("FASTQ files may contain binary information or are compressed")

        #Loop until end of file        
    while rline:

        if rline[0] != '@':
            print(rline)
            raise ValueError("Records in FASTQ files should start with a '@' character. Files may be malformed or out of synch.")

        #Store the existing title, remove right spaces
        title_rline = rline[1:].rstrip()

        #Read another line and remove right spaces, store it (sequence line)
        read1_seq_string = read1_readline().rstrip()

                
        while True:
                        #Read another line, the third of a read entry
            rline = read1_readline()

            if not rline:
                raise ValueError("End of file without quality information. Files may be malformed or out of synch")
            #If the third line starts with a plus, we are finished with the sequence reading
            if rline[0] == '+':
                break
                #Store the sequence without spaces, allow for more than one line of
                        #sequence. Break if a plus is found
            read1_seq_string += rline.rstrip()

        #Store the quality string without right spaces
        read1_quality_string = read1_readline().rstrip()


        while True:
            #Read another line and check whether we have met EOF
            rline = read1_readline()

            if not rline:
                break  # end of file

            if rline[0] == '@' and rline.isalpha() is not True:
                break


            #If we have multiple quality strings, add the next ones, too
            read1_quality_string += rline.rstrip()
            #Return a tuple of title, sequence and quality
        yield (title_rline,  read1_seq_string, read1_quality_string)

        #When we have read the last line, raise built-in exception StopIteration.
        #I don't think it's strictly necessary
    raise StopIteration

# TAG_EXTRACT_FXN
# A function that extracts the tag from a sequence variabale

def tag_extract_fxn(read_seq, blen):
    # This is the function that extracts the UID tags from both the
    # forward and reverse read.  Assigns read1 the sequence from some
    # position to the end, then read2 from some position to the end,
    # then assigns tag1 from the 5'-end to length of the UID tag for
    # read1 and then read 2.
    return(read_seq[:blen])

# HDR_RENAME_FXN
# A function that adds the tag to the title string

def hdr_rename_fxn(read_title, read1_tag):
    # This function renames the header with the formatting of
    # Previous title, *tag from read1*, *tag from read2*,
    # *read designation from original header (for paired reads)*

    #Return the combination of previous title and tag
    return(read_title + ":" + read1_tag)
    




###############################
# MAIN FUNCTION
###############################


def main():

    # Define an argumentparser that interprets command line arguments
    parser =  ArgumentParser()
    parser.add_argument('--infile', dest='infile', help='Path to FASTQ input file. For stdin, use -.', required=True)

    parser.add_argument('--outprefix', dest='outfile', help='Prefix for output files.', required=True)
    parser.add_argument('--taglen', dest='taglen', type=int, default=12,
                        help='Length in bases of the duplex tag sequence.[12]')
    parser.add_argument('--spacerlen', dest='spclen', type=int, default=5,
                        help='Length in bases of the spacer sequence between duplex tag and the start of target DNA. [5]')
    parser.add_argument('--readout', dest='readout', type=int, default=1000000,
                    help='How many reads are processed before progress is reported. [1000000')
    parser.add_argument('--filtspacer', dest='spacer_seq', type=str, default=None,
                        help='Optional: Filter out sequences lacking the inputed spacer sequence. \
                        Not recommended due to significant base calling issues with the invariant spacer sequence')
    parser.add_argument('--tagstats', dest='tagstats', action="store_true",
                        help='Optional: Output tagstats file and make distribution plot of tag family sizes.  \
                        Requires matplotlib to be installed.')

    parser.add_argument('--endtrimming', dest='endtrimming', type=int, default=0, help='Optional: Trim n bases off both ends of the reads [0]')

    
    o = parser.parse_args()


        #Open fastq input file for reading
    if o.infile == '-':
        read1_fastq = sys.stdin
    else:
        read1_fastq = open(o.infile, 'r')

    read1_output = sys.stdout
        
        #Open log file for writing

    logfile = open(o.outfile + '.log',"w") 
    badsequencefile = open(o.outfile + '.badsequences.fq',"w")
    
        #Various counter variables
    readctr = 0    #Read counter
    badspacer = 0
    badspacerno = 0
    goodreads = 0
    badtag = 0
    oldBad = 0
    
    barcode_dict = defaultdict(lambda:0) #Dictionary to store the barcodes
    badspacer_dict = defaultdict(lambda:0) #Dictionary to store bad spacers
    badsequence = False # to distinguis between OK and poor sequences
    #List to store the readlengths in
    readlength = []
    # Generate a SequenceMatcher object, that will be used to find the spacer
    matcher = SequenceMatcher(None,None,o.spacer_seq)
    spclen = o.spclen  # Local variabels that can be changed, if a spacer deviates
    


        
    # FASTQ FILE PROCESSING START

    #Call the generator that iterates through the fastq infile extracting title,
    #sequence and qual
    for read1_title, read1_seq,  read1_qual in fastq_general_iterator(read1_fastq):
        readctr += 1

        #Check whether the spacer is 100% correct. If not, we will accept mistakes except in the three last bases
        if o.spacer_seq != None and (read1_seq[o.taglen:o.taglen + o.spclen] != o.spacer_seq):
            testseq = read1_seq[o.taglen:(o.taglen+o.spclen+3)]  #fetch spacer
            matcher = SequenceMatcher(None,testseq,o.spacer_seq)
            #matcher.set_seq1(testseq)
            blocks = matcher.get_matching_blocks()
            spclen = blocks[-2][0]+blocks[-2][2]  #The last position of the last spacer
                                                #match determines the end of the spacer 
            badspacer=read1_seq[o.taglen:(o.taglen+spclen)] 
                                              
            # If the spacer lacks the three last bases, discard it   
            if badspacer[-3:] != o.spacer_seq[-3:]:
                badsequence = True
                badspacerno += 1
                #Record bad spacers if we are going to report
                if o.tagstats:
                    badspacer_dict[badspacer] += 1

            #End of the if not 100% correct test
            
            #Get the tag from the current sequence read    
        tag1 = tag_extract_fxn(read1_seq, o.taglen)

            #If the tag has only characters and no Ns, add tag to title line
        if (tag1.isalpha() and tag1.count('N') == 0):
            renamed_read1_title =  hdr_rename_fxn(read1_title, tag1)

                #Determine new read sequence
            newread_seq = read1_seq[o.taglen+spclen:]
            newread_qual = read1_qual[o.taglen+spclen:]
            if o.endtrimming > 0:
                newread_seq = newread_seq[o.endtrimming:-o.endtrimming]
                newread_qual = newread_qual[o.endtrimming:-o.endtrimming]

                
                # Write the current read to output file, if its OK
            if not badsequence:
                read1_output.write('@%s\n%s\n+\n%s\n' % (renamed_read1_title, newread_seq, newread_qual))
                #Count successful read converted reads
                goodreads += 1
                    #Increase the counters of the tag dictionary if we will report
                if o.tagstats:
                    barcode_dict[tag1] += 1

                #Register read length
                readlength.append(len(newread_seq))

                #Write the sequence to badsequence file if it has a bad spacer
            elif o.tagstats:
                badsequencefile.write('@%s\n%s\n+\n%s\n' % (read1_title, read1_seq, read1_qual))
                

            
            
            #If the tag was not OK (NNN), count it among the bad tags        
        else:
            badtag += 1

             #Report progress when appropriate
        if readctr % o.readout is 0:
            sys.stderr.write("[TagToHeader] Total sequences processed: %s\n" % readctr)

            # Write error message if all tags were bad in this interval (readout)  
        if badtag == oldBad + o.readout:
            sys.stderr.write("Warning! Potential file error between lines %s and %s." % ((readctr - o.readout) * 4, readctr * 4))
            oldBad = badtag
        badsequence = False   # Re-initalize before new round
        spclen=o.spclen
                
    read1_fastq.close()
    read1_output.close()
    badsequencefile.close()

# PRINT SUMMARY TO LOGFILE
    
    logfile.write("Total sequences processed: %s\n" % readctr)
    logfile.write("Sequences with passing tags: %s\n" % goodreads)
    logfile.write("Bad spacers: %s\n" % badspacerno)
    logfile.write("Bad tags: %s\n" % badtag)
    logfile.write("Mean read length after processing: %d\n" % (sum(readlength)/len(readlength)))

# DO TAG STATISTICS AND PLOTTING
    
    if o.tagstats:
        #Open file for tagstatistics reporting
        tagstatfile = open(o.outfile + '.tagstats', 'w')

        
        #Write the tags and their frequency = number of occurrences
        tagstatfile.write("Tag\tNumber of reads\n")
        for tag in sorted(barcode_dict.keys(),key=barcode_dict.__getitem__,reverse=True):
            tagstatfile.write("%s\t%d\n" % (tag,barcode_dict[tag]))

        #Write the bad spacers and their frequency to file if asked for
        if o.spacer_seq != None:
             badspacerfile = open(o.outfile + '.badspacers', 'w')
             badspacerfile.write("%s is the correct spacer sequence\n\n" % o.spacer_seq)
             badspacerfile.write("Badspacer\tNumber of occurrences\n")
             for badspacer in sorted(badspacer_dict.keys(),key=badspacer_dict.__getitem__,reverse=True):
                 badspacerfile.write("%s\t%d\n" % (badspacer,badspacer_dict[badspacer]))
             badspacerfile.close()

        tagstatfile.close()


        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            familymembers = list(barcode_dict.values())

            #Plot a histogram of tag frequencies
            ax = plt.gca()
            ax.ticklabel_format(style='sci', axis='x')
            plt.hist(familymembers,histtype='stepfilled')

            
            plt.xlabel('Tag family size (reads)')
            plt.ylabel('Frequency')
            plt.savefig(o.outfile + '.pdf')
                        
        except ImportError:
            logfile.write('matplotlib not present. Only tagstats file will be generated.')



if __name__ == "__main__":
    main()
