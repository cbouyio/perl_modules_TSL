package Blast2GFF;

#######################################
#
# A class for tarnsforming blast+ tabular output lines to GFF
#
#######################################

# Currently works with a BLAST+ default result object

use strict;
use warnings;

#use Data::Dumper;

#use vars qw();
#use Exporter 'import';
#our @EXPORT = qw();
#@EXPORT_OK = qw();

#$VERSION = "0.1_alpha";

# Set a global intron length variable! 
our $maxIntronLength = 10000;

sub new {
  # The "constructor"
  my $class = shift;
  my $self = {};
  my $blastResult = shift || die "No BLAST+ result passed\n";
  # If no hit is found return undef so no object is constructed.
  my $blastHit = $blastResult->next_hit || return undef;
  my $blastHSP = $blastHit->next_hsp;
  # Set the BLAST result member variables.
  $$self{'_qseqid'}   = $blastResult->query_name;
  $$self{'_qlength'}  = $blastResult->query_length;
  $$self{'_sseqid'}   = $blastHit->name;
  $$self{'_bitscore'} = $blastHit->bits;
  my $hspNumber  = 1;
  my @hspOrder   = ( $hspNumber );
  my @pIdents    = ( $blastHSP->percent_identity );
  my @hspLengths = ( $blastHSP->length('total') );
  my @qStarts    = ( $blastHSP->start('query') );
  my @qEnds      = ( $blastHSP->end('query') );
  my @sStarts    = ( $blastHSP->start('sbjct') );
  my @sEnds      = ( $blastHSP->end('sbjct') );
  my @eValues    = ( $blastHSP->evalue );
  my @qFrames    = ( $blastHSP->strand('query') );
  my @sFrames    = ( $blastHSP->strand('sbjct') );
  # Swap the query start and end positions in case of reverse hit as GFF3 specs
  # require start to be strictly smaller that the end.
  if ( $blastHSP->start('query') > $blastHSP->end('query')  ) {
    @qStarts = ($blastHSP->end('query'));
    @qEnds   = ($blastHSP->start('query'));
  }
  if ( $blastHSP->start('sbjct') > $blastHSP->end('sbjct')  ) {
    @sStarts = ($blastHSP->end('sbjct'));
    @sEnds   = ($blastHSP->start('sbjct'));
  }
  while ( $blastHSP = $blastHit->next_hsp ) {
    # Check if the current HSP overlaps with any existing one. Procceed only if
    # there is no overlap.
    unless ( hsp_overlaps($blastHSP, \@sStarts, \@sEnds) ) {
      $hspNumber++;
      push @hspOrder, $hspNumber;
      # Swap the query start and end positions in case of reverse hit as GFF3 specs 
      # require start to be strictly smaller that the end.
      if ( $blastHSP->start('sbjct') > $blastHSP->end('sbjct') ) {
        push @sStarts, $blastHSP->end('sbjct');
        push @sEnds, $blastHSP->start('sbjct');
      } else {
        push @sStarts, $blastHSP->start('sbjct');
        push @sEnds, $blastHSP->end('sbjct');
      }
      # Do the same for the query sequence too.
      if ( $blastHSP->start('query') > $blastHSP->end('query') ) {
        push @qStarts, $blastHSP->end('query');
        push @qEnds, $blastHSP->start('query');
      } else {
        push @qStarts, $blastHSP->start('query');
        push @qEnds, $blastHSP->end('query');
      }
      push @pIdents, $blastHSP->percent_identity;
      push @hspLengths, $blastHSP->length('total');
      push @eValues, $blastHSP->evalue;
      push @qFrames, $blastHSP->strand('query');
      push @sFrames, $blastHSP->strand('sbjct');
    }
  }
  # Now sort all the HSPs according to the subject start hit position!
  my @sortedIndices;
  my @sortedSStarts = sort { $a <=> $b } @sStarts;
  for my $start ( @sortedSStarts ) {
    my @sortIndex = grep { $sStarts[$_] eq $start } 0..$#sStarts;
    if ( @sortIndex == 1 ) {
      push @sortedIndices, @sortIndex;
    } else {
      die "Overlaping between HSPs, can not order";
    }
  }
  my @sortedHSPList    = @hspOrder[ @sortedIndices ];
  my @sortedPIdents    = @pIdents[ @sortedIndices ];
  my @sortedHSPLengths = @hspLengths[ @sortedIndices ];
  my @sortedQStarts    = @qStarts[ @sortedIndices ];
  my @sortedQEnds      = @qEnds[ @sortedIndices ];
  my @sortedSEnds      = @sEnds[ @sortedIndices ];
  my @sortedEValues    = @eValues[ @sortedIndices ];
  my @sortedQFrames    = @qFrames[ @sortedIndices ];
  my @sortedSFrames    = @sFrames[ @sortedIndices ];
  # Set the rest BLAST result (HSP related) member variables
  $$self{'_hspNumber'} = $hspNumber;
  $$self{'_hspOrder'}  = \@sortedHSPList;
  $$self{'_pident'}    = \@sortedPIdents;
  $$self{'_hsplength'} = \@sortedHSPLengths;
  $$self{'_qstart'}    = \@sortedQStarts;
  $$self{'_qend'}      = \@sortedQEnds;
  $$self{'_sstart'}    = \@sortedSStarts;
  $$self{'_send'}      = \@sortedSEnds;
  $$self{'_evalue'}    = \@sortedEValues;
  $$self{'_qstrand'}   = \@sortedQFrames;
  $$self{'_sstrand'}   = \@sortedSFrames;

  bless $self, $class;
}


sub hsp_overlaps {
  # Return True if the current HSP overlaps with an existing one
  die "Specify the 3 mandatory arguments" if ( @_ != 3 );
  my $currentHSP = shift;
  my $start = shift;
  my $end = shift;
  my @starts = @$start;
  my @ends = @$end;
  my $currentStart;
  my $currentEnd;
  if ( $currentHSP->start('sbjct') > $currentHSP->end('sbjct') ) {
    $currentStart = $currentHSP->end('sbjct');
    $currentEnd = $currentHSP->start('sbjct');
  } else {
    $currentStart = $currentHSP->start('sbjct');
    $currentEnd = $currentHSP->end('sbjct');
  }
  my $len = @starts;
  for ( my $i = 0; $i < $len; $i++ ) {
    # Condition for overlaping
    if ( ( $currentStart >= $starts[$i] && $currentStart <= $ends[$i] ) || ( $currentEnd >= $starts[$i] && $currentEnd <= $ends[$i] ) ) {
      return 1;
    }
  }
  return 0;
}


sub toGFF {
  # Return the corresponding GFF line.
  my $self = shift;
  die "Specify the 3 mandatory arguments" if ( @_ != 3 );
  my $BLASTtype = $_[0]; # The BLAST program that was used to generate the alignments.
  my $GFFtype   = $_[1]; # The SO (sequence ontology) term of the feature.
  my $coreID        = $_[2]; # The GFF ID of the feature
  my $segments  = "b";
  my $gffLines;
  my $matchID;
  my $segmentID;
  my $flag = 0;
  # Generate an individual GFF3 compatible line for each HSP.
  for (  my $i = 0; $i < $self->hspNumber; $i++ ) {
    # Condition for including the HSP to the same segment
    if ( ( $i > 0 ) and ( ($self->sstart($i) - $self->send($i - 1)) >= $maxIntronLength ) ) {
      $matchID = $coreID . "_" . $segments;
      $segments++;
      $segmentID = $matchID;
      $flag = 1;
    } else {
      if ( $flag ) {
        $matchID = $segmentID;
      } else {
        $matchID = $coreID;
      }
    }
    my $algnsign = $self->get_alignment_sign($i);
    my $gffLine = $self->sseqid . "\t" . $BLASTtype . "\t" . $GFFtype . "\t" . $self->sstart($i) . "\t" . $self->send($i) . "\t" . $self->evalue($i) . "\t" . $algnsign . "\t" . "." . "\t" . "ID=" . $matchID . ";Target=" . $self->qseqid . " " . $self->qstart($i) . " " . $self->qend($i) . "\n";
    $gffLines = $gffLines . $gffLine;
  }
  return $gffLines;
}


sub get_alignment_sign {
  # Return the sign of the match for the GFF file
  my $blastGFF = shift;
  my $i = shift;
  my $algnSign;
  if ( $blastGFF->qstrand($i) * $blastGFF->sstrand($i) > 0 ) {
    $algnSign = "+";
  } elsif ( $blastGFF->qstrand($i) * $blastGFF->sstrand($i) < 0 ) {
    $algnSign = "-";
  } else {
    die "Alignment sign can not have a zero value";
  }
  return $algnSign;
}


# Accessors
sub qseqid {
  my $self = shift;
  return $$self{'_qseqid'};
}

sub qlength {
  my $self = shift;
  return $$self{'_qlength'};
}

sub sseqid {
  my $self = shift;
  return $$self{'_sseqid'};
}

sub bitscore {
  my $self = shift;
  return $$self{'_bitscore'};
}

sub hspNumber {
  my $self = shift;
  return $$self{'_hspNumber'};
}

sub hspOrder {
  my $self = shift;
  my $tmp = $$self{'_hspOrder'};
  my @hspOrder = @$tmp;
  if ( @_ ) {
    return $hspOrder[$_[0]];
  } else {
    return @hspOrder;
  }
}

sub pident {
  my $self = shift;
  my $tmp = $$self{'_pident'};
  my @pident = @$tmp;
  if ( @_ ) {
    return $pident[$_[0]];
  } else {
    return @pident;
  }
}

sub hsplength {
  my $self = shift;
  my $tmp = $$self{'_hsplength'};
  my @hsplength = @$tmp;
  if ( @_ ) {
    return $hsplength[$_[0]];
  } else {
    return @hsplength;
  }
}

#sub gapopen {  my $self = shift;
#  my @args = @_;
#  if ( @args ) {
#    return $$self{'_gapopen'}[$args[0]];
#  } else {
#    return $$self{'_gapopen'};
#  }
#}

sub qstart {
  my $self = shift;
  my $tmp = $$self{'_qstart'};
  my @qstart = @$tmp;
  if ( @_ ) {
    return $qstart[$_[0]];
  } else {
    return @qstart;
  }
}

sub qend {
  my $self = shift;
  my $tmp = $$self{'_qend'};
  my @qend = @$tmp;
  if ( @_ ) {
    return $qend[$_[0]];
  } else {
    return @qend;
  }
}

sub sstart {
  my $self = shift;
  my $tmp = $$self{'_sstart'};
  my @sstart = @$tmp;
  if ( @_ ) {
    return $sstart[$_[0]];
  } else {
    return @sstart;
  }
}

sub send {
  my $self = shift;
  my $tmp = $$self{'_send'};
  my @send = @$tmp;
  if ( @_ ) {
    return $send[$_[0]];
  } else {
    return @send;
  }
}

sub evalue {
  my $self = shift;
  my $tmp = $$self{'_evalue'};
  my @evalue = @$tmp;
  if ( @_ ) {
    return $evalue[$_[0]];
  } else {
    return @evalue;
  }
}

sub qstrand {
  my $self = shift;
  my $tmp = $$self{'_qstrand'};
  my @qstrand = @$tmp;
  if ( @_ ) {
    return $qstrand[$_[0]];
  } else {
    return @qstrand;
  }
}

sub sstrand {
  my $self = shift;
  my $tmp = $$self{'_sstrand'};
  my @sstrand = @$tmp;
  if ( @_ ) {
    return $sstrand[$_[0]];
  } else {
    return @sstrand;
  }
}


# TODO: Posible?? mutators

sub swap_seq_ids {
  # Swap the query ID with the subject ID
  my $self = shift;
  my $tmp = $$self->qseqid;
  $$self{'_qseqid'} = $$self{'_sseqid'};
  $$self{'_sseqid'} = $tmp;
}

sub swap_qstart_qend {
  # Swap the query starting / ending location
  my $self = shift;
  my $tmp = $$self->qstart;
  $$self{'_qstart'} = $$self{'_qend'};
  $$self{'_qend'} = $tmp;
}

sub swap_sstart_send {
  # Swap the subject starting / ending location
  my $self = shift;
  my $tmp = $$self->sstart;
  $$self{'_sstart'} = $$self{'_send'};
  $$self{'_send'} = $tmp;
}

1;


## perldoc stuff...

=head1 NAME

blastPlus2GFF - Module that generates objects representing a GFF3 record out of
NCBI's BLAST+ output (use the default BLAST+ output to get all the supported
features).

=head1 AUTHOR

Costas Bouyioukos costas.bouyioukos@tsl.ac.uk

=head1 SYNOPSIS

    use blastPlus2GFF;
    my $gffLine = new blastPlus2GFF($blastOutLine);

=head1 DESCRIPTION

The parsing object holds some additional (not needed by the GFF) variables...
inherited all from the BLAST parser of the Bio::SeqIO module.

=head1 REQUIRES

No requirements but any user of this package should use the BioPerl BLAST parser
from the Bio::SearchIO package to provide the BLAST result object.

=head1 METHODS

=over

=item new(BioPerl_BLAST_Parsing_Result)

Cary out all the nesseccary parsing and generate a record object holding
information for all the non overlaping HSP.

Return undef if no hits are found in the BLAST output.

=item toGFF(BLAST_program, GFF_type, match_ID)

Return a GFF3 compatible line(s) (as long as the GFF_type argument is a valid SO
name). If multiple HSP have been found a GFF line per individual HSP is
returned. The match ID of this line is kept identical so that all the
information can be represented as an alignment feature in a genome browser. If
the $maxIntornLength (default 5000) has been exceeded the matchID is appended
with an acsending number such that HSPs that are separated will appear as
different alignments in the genome browser.

=back

