package Chemistry::Mok;

$VERSION = '0.22';
# $Id$

use strict;
use warnings;
use Chemistry::Mol;
use Chemistry::File ':auto';
use Chemistry::Pattern;
use Chemistry::Bond::Find qw(find_bonds assign_bond_orders);
use Text::Balanced ':ALL';
use Scalar::Util 'blessed';

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
            { 'Chemistry::Mok::Opts' => 
                qr/[gopGOP]+/ },
        ],
    );
    die "Mok: error extracting: $@" if $@;
    #use Data::Dumper; print Dumper \@toks; print "$code\n";
    @toks;
}

sub parse {
    my (@toks)  = @_;

    my (@subs, @blocks);
    for my $tok (@toks) {
        blessed $tok or die "unparsable token '$tok'\n";
    }
    while (my $tok = shift @toks) {
        #print "tok = $$tok\n";
        if ($tok->isa("Chemistry::Mok::Sub")) {
            my $next = shift @toks 
                or die "unexpected end of mok program after $$tok\n";
            if ($next->isa("Chemistry::Mok::Block")) {
                push @subs, "$$tok $$next";
            } else {
                die "unexpected token $$tok after Sub; expected Block\n";
            }
        } elsif ($tok->isa("Chemistry::Mok::Patt")) {
            my $next = shift @toks
                or die "unexpected end of mok program after $$tok\n";
            my $opts = '';
            if ($next->isa("Chemistry::Mok::Opts")) {
                $opts = $$next;
                $next = shift @toks
                    or die "unexpected end of mok program after $$tok\n";
            }
            if ($next->isa("Chemistry::Mok::Block")) {
                push @blocks, { patt => $$tok, opts => $opts, block => $$next};
            } else {
                die "unexpected token $$tok after Patt; expected Block\n";
            }
        } elsif ($tok->isa("Chemistry::Mok::Block")){
            push @blocks, { patt => '', opts => '', block => $$tok};
        } elsif ($tok->isa("Chemistry::Mok::Comment")){
            # do nothing
        } else {
            die "unexpected token $$tok\n";
        }
    }
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
            $patt = Chemistry::Pattern->parse($patt_str, format => $format);
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

The name of the format which will be used for parsing slash-delimited patterns.
Mok versions until 0.16 only used the 'smiles' format, but newer versions can
use the 'smarts' format as well.

=back

=cut

sub new {
    my ($class, $code, @a) = @_;
    my %opts;
    unshift @a, "package" if (@a == 1);
    %opts = @a;
        
    my $usr_pack = $opts{package} || "Default"; 
    my $format   = $opts{pattern_format} || "smiles"; 

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

=item mol_class

The molecule class used for reading the files. Defaults to Chemistry::Mol.

=item format

The format used when calling $mol_class->read. If not given, $mol_class->read
tries to identify the format automatically.

=item find_bonds

If set to a true value, find bonds. Use it when reading files with no bond
information but 3D coordinates to detect the bonds if needed (for example, if
you want to do match a pattern that includes bonds). If the file has explicit
bonds, mok will not try to find the bonds, but it will reassign the bond orders
from scratch.

=back

=cut

sub run {
    my ($self, $opt, @args) = @_;
    # MAIN LOOP
    my $mol_class = $opt->{mol_class} || "Chemistry::Mol";
    FILE: for my $file (@args) {
        my (@mols) = $mol_class->read($file, format => $opt->{format},
            mol_class => $opt->{mol_class});
        MOL: for my $mol (@mols) {
            if ($opt->{find_bonds}) {
                find_bonds($mol) unless $mol->bonds;
                assign_bond_orders($mol);
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


__END__

=back

=head1 VERSION

0.22

=head1 SEE ALSO

L<mok>, L<http://www.perlmol.org/>

=head1 AUTHOR

Ivan Tubert-Brohman E<lt>itub@cpan.orgE<gt>

=head1 COPYRIGHT

Copyright (c) 2004 Ivan Tubert-Brohman. All rights reserved. This program is
free software; you can redistribute it and/or modify it under the same terms as
Perl itself.

=cut

