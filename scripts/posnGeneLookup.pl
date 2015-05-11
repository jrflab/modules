#!/usr/bin/perl
# lookup the gene(s) at a base position (first column of the input)

use strict;

use DBI;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;

print STDERR "Connecting to Ensembl core...\n";
my $dbCore = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => "anonymous",
    -dbname => "homo_sapiens_core_65_37",
    -host   => "localhost",
    -pass   => "",
    -driver => 'mysql',
    -port   => 33387
);
print STDERR "Connected.\n";

my $slice_adaptor = $dbCore->get_SliceAdaptor();
while (<>) {
    chomp;
    my @F = split / /;
    print STDERR "Looking up position $F[0]\n";
    my ($chr, $pos) = split /:/, $F[0];
    $chr =~ s/chr//;
    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $pos, $pos);
    my @genes = @{$slice->get_all_Genes()};
    my @ids;
    my @strands;
    while (my $gene = shift @genes) {
        my $stable_id = $gene->stable_id();
        my $strand = $gene->strand();
        push @ids, $stable_id;
        push @strands, $strand;
    }
    print "$chr:$pos " . ((@ids > 0)? join("|", @ids) . ' ' . join("|", @strands): "NA NA") . ' ' . join(" ", @F[1..$#F]) . "\n";
}
