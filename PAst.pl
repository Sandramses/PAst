#!/usr/bin/perl -w
use strict;

#the Pseudomonas aeruginosa serotyper (PAst) v3
    
if ($ARGV[0] eq "help" or $ARGV[0] eq "-h" or $ARGV[0] eq "man") {
    print   "\nProgram:\tPAst (Pseudomonas aeruginosa serotyper)\nVersion:\t1.0 (using NCBI blast+ 2.2.29)\n\nUsage:\t\t./PAst.pl <path/to/BLASTbin> <path/to/output/directory> <path/to/input/directory> <path/to/OSAdatabase>\n\nBLAST bin:\tPath to the bin containing blasts executables (ex. /usr/bin/ncbi-blast-2.2.29+/bin/)\nOutput dir:\tDirectory where blast reports, OSA multifasta files and result summary (serotyping.txt) are placed\nInput dir:\tDirectory containing assembled input genomes in multifasta format (and nothing else)\nOSA database:\tPath to OSA database file (downloadable from Github)\n\n";
} else {
    #Configuration (user specified):
    my $BLASTbin = $ARGV[0];
    my $outdir = $ARGV[1];
    my $inputdir = $ARGV[2];
    my $db = $ARGV[3];
    
    #####################################################
    ################# Input acquisition #################
    #####################################################
    my @IDs = [];
    
    #input filenames are retrieved from $inputdir
    opendir DIR, $inputdir or die "Cannot open input directory $inputdir: $!\n";
    @IDs = readdir DIR;
    closedir DIR;
    
    #directory references are removed from list
    my @sorted = sort @IDs;
    @IDs = @sorted;
    splice(@IDs, 0, 2);
    
    #number of files localized is printed for user and pipeline initiated
    print "There are ", scalar(@IDs), " files in the input directory, these are now processed by PAST.\n";
    
    #####################################################
    ################# Serotype pipeline #################
    #####################################################
    
    my @output = ();
    
    #each file in the input directory is processed and typed separately
    foreach my $file (@IDs) {
        
        #blastn analysis is performed in terminal, output saved in $outdir
        system("cd $outdir");
        system("$BLASTbin/blastn -query $inputdir/$file -subject $db -out $outdir/BLAST$file");
        
        #information about OSA clusters (sequence and length) is extracted from blastn output
        my %data = &extractBLAST_sub($file);
        
        #the isolate is typed if data have been extracted from blast output, otherwise error messages are printed
        if (exists $data{O1}) {
            #in silico serotype is determined based on blastn data
            my @serotype = &serotype_sub(\%data);
            
            #result summary for file is saved in @output and individual file output is generated (sub)
            push(@output, &output_sub(\$file, \%data, \@serotype));
        } else {
            #error messages if data was not extracted from blastn output
            print "- PAst encountered problems with BLAST analysis of $file, please check manually.\n";
            push(@output, "No data could be extracted from the BLAST output for $file, please check file manually.")
        }
    }
    
    #####################################################
    ################# Output generation #################
    #####################################################
    
    #combined result summary for all input files is generated as serotyping.txt
    open(OUT, '>', "$outdir/serotyping.txt") or die "Output file serotyping.txt could not be written: $!";
    
    #header line is followed by result summaries for individual files (tab-delimited)
    print OUT "filename\tserotype\t%coverage of OSA\tOSA fragments\n";
    print OUT join("\n", @output);
        
    close OUT;
    
    #user is informed the analysis is finalized an directed to the output files
    print "The analysis has finished. Results (serotyping.txt) and BLAST reports are placed in $outdir\n";
    
    #####################################################
    #################### Subroutines ####################
    #####################################################
    
    #subroutine for extracting data from blastn output
    sub extractBLAST_sub {
        #open blast output file
        open(IN, '<', "$outdir/BLAST$_[0]") or print "ERROR blastn output file BLAST$_[0] could not be opened or was not generated: $!";
        
        #variables
        my %extract = ();
        my ($Otype, $cov) = '';
        
        #parse blastn output and extract informatin into %extract
        while (defined(my $line = <IN>)) {
            #extract names of OSA clusters with hits in query        
            $Otype = $1 if $line =~ m/^Subject=\s+([^\s]+)/;
            
            #extract %identity of blast hit
            $cov = $1 if $line =~ m/^\s+Identities\s+=\s+\d+[^\d]\d+[^\d]+(\d+)/;
            
            #extract dna sequence and length if %identity exceeds 94%
            if ($line =~ m/^Query\s+(\d+)\s+(\w+)\s+(\d+)/) {    
                $extract{$Otype}{sequence} .= uc($2) unless $cov < 95;
                $extract{$Otype}{length} += length($2) unless $cov < 95;
            }
                
            #add ' ' after each type, to separate multiple hits
            $extract{$Otype}{sequence} .= ' ' if $line =~ m/^\s+Score/;
            
        }    
        close IN;
        
        #return extracted data
        return %extract;
    }
    
    #subroutine for determining in silico serotype    
    sub serotype_sub {
        #variables
        my @type = ();
        my $besttype = '';
        my $tmp = 0;
        my $flag = 0;
        
        #OSA ref sizes
        my %OSAref = (O1 => 18368, O2 => 23303, O3 => 20210, O4 => 15279, O6 => 15649, O7 => 19617, O9 => 17263, O10 => 17635,
                      O11 => 13868, O12 => 25864, O13 => 14316, WzyB => 1140);
        
        #detection of WzyB
        foreach my $O (sort keys %{$_[0]}) {
            if ($O eq "WzyB") {
                #calculate coverage of wzyB
                ${$_[0]}{$O}{coverage} = 100 - (($OSAref{$O} - ${$_[0]}{$O}{length})/$OSAref{$O} * 100);
                
                #flag is on if wzyB was detected
                $flag = 1 if ${$_[0]}{$O}{coverage} > 95;
            }
        }
        
        #assignment of in silico serotype
        foreach my $O (sort keys %{$_[0]}) {
            
            #calculate coverage of OSA gene clusters
            ${$_[0]}{$O}{coverage} = 100 - (($OSAref{$O} - ${$_[0]}{$O}{length})/$OSAref{$O} * 100);
            
            #save serotypes with over 95% coverage in @type and best hit
            if ($O ne "WzyB") {
                if (${$_[0]}{$O}{coverage} > 95) {
                
                    #O2 + wzyB = O2/O16
                    if ($flag == 1 and $O eq "O2") {
                        push(@type, $O);
                    #O2 - wzyB = O5/O18/O20
                    } elsif ($flag == 0 and $O eq "O2") {
                        push(@type, "O5");
                        ${$_[0]}{O5}{coverage} = ${$_[0]}{$O}{coverage};
                    #other types based only on OSA
                    } else {push(@type, $O);}
                
                #no significant hit, best type is stored
                } elsif ($tmp < ${$_[0]}{$O}{coverage}) {
                    $tmp = ${$_[0]}{$O}{coverage};
                    $besttype = $O;
                }            
            }
        }
        
        #if no significant hit is identified the best is added to @type
        if (!defined $type[0]) {
            push(@type, $besttype);
        }
        
        #return identified serotype of isolates
        return (@type);
    }
    
    #subroutine for generating output
    sub output_sub {
        #variable
        my @results = ();
        
        #for each assigned serotype of the isolate print message to user about analysis result (one type, indecisive or no type)
        foreach my $type (@{$_[2]}) {
            
            #print serotyping result for user
            if (${$_[1]}{$type}{coverage} eq "NT") {
                print "- No serotype could be detected for isolate ${$_[0]}.\n";
            } elsif (scalar(@{$_[2]}) == 1 and ${$_[1]}{$type}{coverage} > 95) {
                print "- The serotype of the isolate ${$_[0]} is ${$_[2]}[0].\n";
            } elsif (scalar(@{$_[2]}) > 1) {
                print "- The serotype of the isolate ${$_[0]} is indecisive.\n";
            }
            
            #generate fasta file with sequence of OSA cluster
            open(OUT, '>', "$outdir/$type-OSA_${$_[0]}") or print "ERROR output fasta file $outdir/$type-OSA_${$_[0]} could not be generated: $!\n";
            
            my @seqs = ();
            #split sequences from different blastn hits
            if ($type eq 'O5') {
                @seqs = split(" ", ${$_[1]}{O2}{sequence});
            } else {
                @seqs = split(" ", ${$_[1]}{$type}{sequence});
            }
            
            my $i = 1;
            #print header and sequence (60 chars per line) for each blastn hot OSA fragment    
            foreach my $seq (@seqs) {
                print OUT "> $type-OSA.$i, size ", length($seq), "\n";
                for (my $j = 0; $j < length($seq); $j += 60) {
                    print OUT substr($seq, $j, 60), "\n";
                }
                $i++;
            }
            
            close OUT;
            
            #generate result summary for each assigned serotype
            push(@results, "${$_[0]}\t$type\t${$_[1]}{$type}{coverage}\t".scalar(@seqs)."\n");
        }
        
        #return result summary
        return @results;
    }
}
__END__
