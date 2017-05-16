
# Check for correct PERL version (>5.0)

if ($] < 5.0) {
    print STDERR "\nFATAL:  You have specified Perl version ",
    $]+0.00," in the Makefile,\n        you must use version 5.0 or greater.\n\n",
    "Install Perl 5.0 or greater and/or reset the \$PERLDIR variable in\n",
    "'Makefile' to the directory in which Perl 5.0 or greater is installed.\n\n";
    exit (1);
}

exit (0);

