#!/usr/bin/perl
use Cwd;



#qsub -N merge -P PR215 -q serial -o OUT/5.out -e ERR/5.err -l select=1:ncpus=1:mem=6G -l walltime=02:00:00 ./4c_make_mapping_stats.pl

my $dir = getcwd;

print $dir;

open FPIN, "<".$dir."/../resources/param.txt" or die $dir."/../resources/param.txt";
while (<FPIN>) {
	if (/^sampleNames\t([^\n]+)\n$/) {
		my $line = $1;
		$line =~ s/^ //g;
		$line =~ s/ $//g;
		@samples = split / /, $line;
	}
}
close FPIN;

## build data structure

my $ds;

foreach my $sample (@samples) {
        open FPIN, "<".$dir."/../output/".$sample."/fastqc/untrimmed/".$sample."_1_fastqc.html" or die;

        while (<FPIN>) {

                if (/<td>Total Sequences<\/td><td>(\d+)<\/td>/) {
                        $ds->{$sample}->{"total"} = ($1 * 2);
                        $ds->{$sample}->{"ptotal"} = "100";
                } 
	}

        close FPIN;
}



foreach my $sample (@samples) {
	open FPIN, "<".$dir."/../output/".$sample."/bamtools/".$sample.".markdup.stat.txt" or die;

	while (<FPIN>) {

		if (/^Total reads:\s+(\d+)\n$/) {
			$ds->{$sample}->{"trimtotal"} = $1;
			$ds->{$sample}->{"ptrimtotal"} = (($1 / $ds->{$sample}->{"total"})*100);
		} elsif (/^Mapped reads:\s+(\d+)[^\(]+\(([^\%]+)[\%]/) {
                        $ds->{$sample}->{"mapped"} = $1;
                        $ds->{$sample}->{"pmapped"} = $2;
                } elsif (/^Forward strand:\s+(\d+)[^\(]+\(([^\%]+)[\%]/) {
                        $ds->{$sample}->{"forward"} = $1;
                        $ds->{$sample}->{"pforward"} = $2;
                } elsif (/^Reverse strand:\s+(\d+)[^\(]+\(([^\%]+)[\%]/) {
                        $ds->{$sample}->{"reverse"} = $1;
                        $ds->{$sample}->{"preverse"} = $2;
                } elsif (/^Duplicates:\s+(\d+)[^\(]+\(([^\%]+)[\%]/) {
                        $ds->{$sample}->{"duplicate"} = $1;
                        $ds->{$sample}->{"pduplicate"} = $2;
                }
	}

	close FPIN;

	my $wc = `wc -l /scratch/mtera/analysis/an0099/bin/../output/$sample\/bedtools/$sample\.ontarget.sam`;	
	die unless ($wc =~ /^(\d+)\s/);
	$ds->{$sample}->{"target"} = $1;
        $ds->{$sample}->{"ptarget"} = (($1 / $ds->{$sample}->{"mapped"})*50);

}



open FPOUT, ">".$dir."/../output/mapping_stats.txt";

print FPOUT "sampleID\ttotalReadPairs\ttotalReadPairs (\%)\ttotalReadPairsPostTrim\ttotalReadPairsPostTrim (\%)\tmappedReadPairs\tmappedPairReads (\%)\t";
print FPOUT "forwardReads\tforwardReads (\% of mapped)\treverseReads\treverseReads (\% of mapped)\t";
print FPOUT "duplicateReadPairs\tduplicateReadPairs (\% of mapped)\tonTargetReadPairs\tonTargetReadPairs (\% of mapped)\n";

foreach my $sample (@samples) {
	print FPOUT $sample;
	print FPOUT "\t".($ds->{$sample}->{"total"} / 2);
	print FPOUT "\t".$ds->{$sample}->{"ptotal"};
	print FPOUT "\t".($ds->{$sample}->{"trimtotal"} / 2);
	print FPOUT "\t".$ds->{$sample}->{"ptrimtotal"};
	print FPOUT "\t".($ds->{$sample}->{"mapped"});
	print FPOUT "\t".$ds->{$sample}->{"pmapped"};
	print FPOUT "\t".$ds->{$sample}->{"forward"};
	print FPOUT "\t".$ds->{$sample}->{"pforward"};
	print FPOUT "\t".$ds->{$sample}->{"reverse"};
	print FPOUT "\t".$ds->{$sample}->{"preverse"};
	print FPOUT "\t".($ds->{$sample}->{"duplicate"});
	print FPOUT "\t".$ds->{$sample}->{"pduplicate"};
	print FPOUT "\t".($ds->{$sample}->{"target"} / 2);
	print FPOUT "\t".$ds->{$sample}->{"ptarget"};
	print FPOUT "\n";

}
