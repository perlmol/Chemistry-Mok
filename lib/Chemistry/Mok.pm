package Chemistry::Mok;

$VERSION = '0.24';
# $Id$

use strict;
use warnings;
use Chemistry::Mol;
use Chemistry::File ':auto';
use Chemistry::Pattern;
use Chemistry::Bond::Find qw(find_bonds assign_bond_orders);
use Chemistry::Ring 'aromatize_mol';
use Text::Balanced ':ALL';
use Scalar::Util 'blessed';
use Data::Dumper;
use Carp;

our $DEBUG = 0;

=head1 NAME

Chemistry::Mok - molecular awk interpreter

=head1 SYNOPSIS

    use Chemistry::Mok;
    $code = '/CS/g{ $n++; $l += $match->bond_map(0)->length }
        END { printf "Average C-S bond length: %.3f\n", $l/$n; }';

    my $mok = Chemistry::Mok->new($code);
    $mok->run({ format => mdlmol }, glob("*.mol"));

=head1 DESCRIPTION

This module is the engine behind the mok program. See mok(1) for a detailed
description of the language. Mok is part of the PerlMol project,
L<http://www.perlmol.org>.

=head1 METHODS

=over

=cut

sub tokenize {
    my ($code) = @_;

    $code =~ s/\s*$//;
    $code =~ s/^\s*#.*//g; # remove comments at the top of the file
    unless($code =~ /^\s*([\/{#]|sub|BEGIN|END)/) {
        $code = "{$code}"; # add implicit brackets for simple one-liners
    }
    #print "code = '$code'\n";
    # (patt opt?)? code | sub code
    my @toks = extract_multiple(my $c = $code,
        [
            { 'Chemistry::Mok::Comment' => 
                qr/\s*#.*\s*/ },
            { 'Chemistry::Mok::Patt' => 
                sub { scalar extract_delimited($_[0],'/') } },
            { 'Chemistry::Mok::Sub'  => 
                qr/\s*(END|BEGIN|sub\s\w+)\s*/ },
            { 'Chemistry::Mok::Block' => 
                sub { scalar extract_codeblock($_[0],'{') } },
            { 'Chemistry::Mok::PattLang' => 
                qr/\s*(\w+):(?=\s*\/)/ },
            { 'Chemistry::Mok::Opts' => 
                qr/[gopGOP]+/ },
        ],
    );
    die "Mok: error extracting: $@" if $@;
    print "TOKENS:\n", Dumper(\@toks), "\nCODE:<<<<$code>>>>\n\n"
        if $DEBUG;
    @toks;
}

sub parse {
    my (@toks)  = @_;

    my (@subs, @blocks);
    for my $tok (@toks) {
        blessed $tok or die "unparsable token '$tok'\n";
    }

###  new parser

    my $st = 1;
    my ($patt, $opts, $block, $sub, $pattlang) = ('') x 5;
    my ($save) = 0;
    while (my $tok = shift @toks) {
        next if $tok->isa("Chemistry::Mok::Comment");
        if ($st == 1) {
            if ($tok->isa("Chemistry::Mok::Block")){
                $block = $$tok,     $save = 1;
            } elsif ($tok->isa("Chemistry::Mok::Sub")) {
                $sub = $$tok,       $st = 5,    next;
            } elsif ($tok->isa("Chemistry::Mok::PattLang")) {
                $pattlang = $$tok,  $st = 4,    next;
            } elsif ($tok->isa("Chemistry::Mok::Patt")) {
                $patt = $$tok,      $st = 2,    next;
            }
        } elsif ($st == 2) {
            if ($tok->isa("Chemistry::Mok::Block")){
                $block = $$tok,     $save = 1;
            } elsif ($tok->isa("Chemistry::Mok::Opts")){
                $opts = $$tok,      $st = 3,    next;
            }
        } elsif ($st == 3) {
            if ($tok->isa("Chemistry::Mok::Block")){
                $block = $$tok,     $save = 1;
            }
        } elsif ($st == 4) {
            if ($tok->isa("Chemistry::Mok::Patt")){
                $patt = $$tok,      $st = 2,    next;
            }
        } elsif ($st == 5) {
            if ($tok->isa("Chemistry::Mok::Block")){
                $block = $$tok,     $save = 1;
            }
        } else {
            confess "unknown state '$st'";
        }
        if ($save) { # save block and go back to state 1
            if ($sub) {
                push @subs, "$sub $$tok";
            } else {
                push @blocks, { patt => $patt, opts => $opts, 
                    pattlang => $pattlang, block => $$tok};
            }
            $patt = $opts = $pattlang = $block = $sub = '';
            $st = 1,    $save = 0,  next;
        } else {
            die "unexpected token '$$tok' (type '" . ref($tok) . "'\n";
        }
    }
    print "BLOCKS\n", Dumper(\@blocks), "\nSUBS:\n", Dumper(\@subs), "\n"
        if $DEBUG;

    \@subs, \@blocks;
}

sub compile_subs {
    my ($pack, @subs) = @_;
    for my $sub (@subs) {
        eval <<END;
            package Chemistry::Mok::UserCode::$pack;
            no strict;
            no warnings;
            $sub
END
        die "Mok: error compiling sub: $@" if $@;
    }
}

sub compile_blocks {
    my ($pack, $format, @blocks) = @_;
    my @compiled_blocks;

    for my $block (@blocks) {
        #use Data::Dumper; print Dumper $block;
        my $code = $block->{block};
        my $sub = eval <<END;
            package Chemistry::Mok::UserCode::$pack;
            no strict;
            no warnings;
            sub {
                my (\$mol, \$file, \$match, \$patt) = \@_;
                my (\$MOL, \$FILE, \$MATCH, \$PATT) = \@_;
                my (\@A) = \$MATCH ? \$MATCH->atom_map : \$MOL->atoms;
                my (\@B) = \$MATCH ? \$MATCH->bond_map : \$MOL->bonds;
                $block->{block};
            }
END
        die "Mol: Error compiling block: $@" if $@;

        my ($patt, $patt_str);
        if ($block->{patt}) {
            $block->{patt} =~ m#^/(.*)/$#;
            $patt_str = $1;
            $patt = Chemistry::Pattern->parse($patt_str, 
                format => $block->{pattlang} || $format);
            $patt->attr(global => 1) if $block->{opts} =~ /g/;
            $patt->options(overlap => 0) if $block->{opts} =~ /O/;
            $patt->options(permute => 1) if $block->{opts} =~ /p/;
        } 
        push @compiled_blocks, {'sub' => $sub, 
            patt => $patt, patt_str => $patt_str};
    }
    \@compiled_blocks;
}

=item Chemistry::Mok->new($code, %options)

Compile the code and return a Chemistry::Mok object. Available options:

=over

=item C<package>

If the C<package> option is given, the code runs in the
Chemistry::Mok::UserCode::$options{package} package instead of the
Chemistry::Mok::UserCode::Default package. Specifying a package name is
recommended if you have more than one mok object and you are using global
varaibles, in order to avoid namespace clashes.

=item C<pattern_format>

The name of the format which will be used for parsing slash-delimited patterns
that don't define an explicit format. Mok versions until 0.16 only used the
'smiles' format, but newer versions can use other formats such as 'smarts',
'midas', 'formula_pattern', and 'sln', if available. The default is 'smarts'.

=back

=cut

sub new {
    my ($class, $code, @a) = @_;
    my %opts;
    unshift @a, "package" if (@a == 1);
    %opts = @a;
        
    my $usr_pack = $opts{package} || "Default"; 
    my $format   = $opts{pattern_format} || "smarts"; 

    # import convenience functions into the user's namespace
    eval <<EVAL;
          package Chemistry::Mok::UserCode::$usr_pack;
          Chemistry::Atom->import(':all');
          Math::VectorReal->import(':all');
          sub println { print "\@_", "\n" }
EVAL
    my @toks = tokenize($code);
    my ($subs, $blocks) = parse(@toks);
    compile_subs($usr_pack, @$subs);
    my $mok = compile_blocks($usr_pack, $format, @$blocks);
    bless $mok, ref $class || $class;
}

=item $mok->run($options, @args)

Run the code on the filenames contained in @args. $options is a hash reference
with runtime options. Available options:

=over

=item aromatize          

"Aromatize" each molecule as it is read. This is needed for example for
matching SMARTS patterns that use aromaticity or ring primitives.

=item delete_dummies

Delete dummy atoms after reading each molecule. A dummy atom is defined as an
atom with an unknown symbol (i.e., it doesn't appear on the periodic table), or
an atomic number of zero.

=item find_bonds

If set to a true value, find bonds. Use it when reading files with no bond
information but 3D coordinates to detect the bonds if needed (for example, if
you want to do match a pattern that includes bonds). If the file has explicit
bonds, mok will not try to find the bonds, but it will reassign the bond orders
from scratch.

=item format

The format used when calling $mol_class->read. If not given, $mol_class->read
tries to identify the format automatically.

=item mol_class

The molecule class used for reading the files. Defaults to Chemistry::Mol.

=back

=cut

sub run {
    my ($self, $opt, @args) = @_;
    # MAIN LOOP
    my $mol_class = $opt->{mol_class} || "Chemistry::Mol";
    FILE: for my $file (@args) {
        my (@mols) = $mol_class->read(
            $file, 
            format      => $opt->{format},
            mol_class   => $opt->{mol_class},
        );
        MOL: for my $mol (@mols) {
            if ($opt->{delete_dummies}) {
                $_->delete for grep { ! $_->Z } $mol->atoms;
            }
            if ($opt->{find_bonds}) {
                find_bonds($mol) unless $mol->bonds;
                assign_bond_orders($mol);
            }
            if ($opt->{aromatize}) {
                aromatize_mol($mol);
            }
            BLOCK: for my $block (@$self) {
                my ($code_block, $patt, $patt_str) = 
                    @{$block}{qw(sub patt patt_str)};
                if ($patt) {
                    MATCH: while ($patt->match($mol)) {
                        $code_block->($mol, $file, $patt, $patt_str);
                        last unless $patt->attr('global');
                    }
                } else {
                    $code_block->($mol, $file, $patt, $patt_str);
                }
            }
        }
    }
}

1;

__END__

=back

=head1 VERSION

0.24

=head1 SEE ALSO

L<mok>, L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2005 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

