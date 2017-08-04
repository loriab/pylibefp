import sys
import math


def success(label):
    """Function to print a '*label*...PASSED' line to screen.
    Used by :py:func:`util.compare_values` family when functions pass.

    """
    msg = '\t{0:.<66}PASSED'.format(label)
    print(msg)
    sys.stdout.flush()


def failed(label):
    print('\t{0:.<66}FAILED'.format(label))


class ComparisonTestError(Exception):
    """Error called when a test case fails due to a failed
    compare_values() call. Prints error message *msg* to standard
    output stream and output file.

    """
    def __init__(self, msg):
        PsiException.__init__(self, msg)
        self.message = '\nPsiException: %s\n\n' % msg


# Test functions
def compare_values(expected, computed, digits, label, exitonfail=True):
    """Function to compare two values. Prints :py:func:`util.success`
    when value *computed* matches value *expected* to number of *digits*
    (or to *digits* itself when *digits* < 1 e.g. digits=0.04). Performs
    a system exit on failure unless *exitonfail* False, in which case
    returns error message. Used in input files in the test suite.

    """
    if digits > 1:
        thresh = 10 ** -digits
        message = ("\t%s: computed value (%.*f) does not match (%.*f) to %d digits." % (label, digits+1, computed, digits+1, expected, digits))
    else:
        thresh = digits
        message = ("\t%s: computed value (%f) does not match (%f) to %f digits." % (label, computed, expected, digits))
    if abs(expected - computed) > thresh:
        print(message)
        if exitonfail:
            return False
            #raise TestComparisonError(message)
    if math.isnan(computed):
        print(message)
        print("\tprobably because the computed value is nan.")
        if exitonfail:
            return False
            #raise TestComparisonError(message)
    success(label)
    return True

