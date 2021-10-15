#!/usr/bin/perl

my @samples;
open FPIN, "<../resources/param.txt" or die;
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

	my $count = 0;

	open FPIN, "<../output/".$sample."/featurecount/".$sample.".markdup.featurecount" or die;
	<FPIN>;
	<FPIN>;
	while (<FPIN>) {
		if (/^([^\t]+)\t.+\t([^\t\n]*)\n$/) {
			my ($id, $val) = ($1, $2);
			die $id if (exists $ds->{$id}->{$sample});
			$ds->{$id}->{$sample} = $val;
			$count += $val;
		}
	}	
	close FPIN;

	print "markdup: ".$sample." ".$count."\n";	
}

open FPOUT, ">../output/all.markdup.featurecount" or die;
print FPOUT "tracking_id";
foreach my $sample (@samples) { print FPOUT "\t".$sample; }
print FPOUT "\n";

foreach my $id (sort keys %{$ds}) {
	print FPOUT $id;

	foreach my $sample (@samples) {	

		print FPOUT "\t".$ds->{$id}->{$sample};
	}
	print FPOUT "\n";
}
close FPOUT;

my $ds;
foreach my $sample (@samples) {

	my $count = 0;

        open FPIN, "<../output/".$sample."/featurecount/".$sample.".rmdup.featurecount" or die;
	<FPIN>;
	<FPIN>;
        while (<FPIN>) {
                if (/^([^\t]+)\t.+\t([^\t\n]*)\n$/) {
                        my ($id, $val) = ($1, $2);
                        die $id if (exists $ds->{$id}->{$sample});
                        $ds->{$id}->{$sample} = $val;
			$count += $val;
                }
        }
	close FPIN;

	print "rmdup: ".$sample." ".$count."\n";  
}

open FPOUT, ">../output/all.rmdup.featurecount" or die;
print FPOUT "tracking_id";
foreach my $sample (@samples) { print FPOUT "\t".$sample; }
print FPOUT "\n";

foreach my $id (sort keys %{$ds}) {
        print FPOUT $id;

        foreach my $sample (@samples) {

                print FPOUT "\t".$ds->{$id}->{$sample};
        }
	print FPOUT "\n";
}
close FPOUT;



