import sys
import math
import pprint
import deepdiff


def _success(label):
    """Function to print a '*label*...PASSED' line to screen.
    Used by :py:func:`util.compare_values` family when functions pass.

    """
    msg = '\t{0:.<66}PASSED'.format(label)
    print(msg)
    sys.stdout.flush()


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
    if math.isnan(computed):
        print(message)
        print("\tprobably because the computed value is nan.")
        if exitonfail:
            return False
    _success(label)
    return True


def compare_integers(expected, computed, label, exitonfail=True):
    """Function to compare two integers. Prints :py:func:`util.success`
    when value *computed* matches value *expected*.
    Performs a system exit on failure. Used in input files in the test suite.

    """
    if (expected != computed):
        print("\t%s: computed value (%d) does not match (%d)." % (label, computed, expected))
        return False
    _success(label)
    return True


def compare_dicts(expected, computed, tol, label):
    """Compares dictionaries *computed* to *expected* using DeepDiff

    Float comparisons made to *tol* significant decimal places.

    Note that a clean DeepDiff returns {}, which evaluates to False, hence the compare_integers.

    """
    ans = deepdiff.DeepDiff(expected, computed, significant_digits=tol, verbose_level=2)
    clean = not bool(ans)
    if not clean:
        pprint.pprint(ans)
    return compare_integers(True, clean, label)

