use 5.006001;
use ExtUtils::MakeMaker;

WriteMakefile(
    'NAME'		=> 'Chemistry::Mok',
    'VERSION_FROM'	=> 'lib/Chemistry/Mok.pm', 
    'PREREQ_PM'		=> {
        'Chemistry::Mol' => 0.1,
        'Chemistry::Pattern' => 0.1,
        'Chemistry::File::SMILES' => 0.21,
        'Test::Simple' => 0,
        'Scalar::Util' => 1.01,
        'Text::Balanced' => 1.87,
        'Math::VectorReal' => 1.0,
    },
    'EXE_FILES' => ['mok'],
);
