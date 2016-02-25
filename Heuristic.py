"""
Module docstring goes here, with an explanation of why one needs to
use Heuristic objects instead of just calling functions.
"""
import abc
import unittest


class Result(object):
    """
    The result of a call of Heuristic.run().  This is an abstract base
    class with only two subclasses: Success and Failure.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def is_success(self):
        """Return True iff this is a Success instance."""
        pass

    def is_failure(self):
        return not self.is_success()


class Success(Result):
    """
    A result of a successful call to Heuristic.run().  Contains the
    result value.
    """
    def __init__(self, value):
        self.value = value

    def __repr__(self):
        return 'Success(%r)' % self.value

    def __str__(self):
        return self.__repr__()

    def is_success(self):
        return True


class Failure(Result):
    """
    A result of a failed call of Heuristic.run().  Contains the
    exceptions raised by the calculation.
    """
    def __init__(self, exceptions):
        self.exceptions = exceptions

    def __repr__(self):
        return 'Failure(%r)' % self.exceptions

    def __str__(self):
        return 'Failure'

    def is_success(self):
        return False


class Heuristic(object):
    """
    A composable representation of a calculation on one argument that
    may fail.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def run(self, arg):
        """
        Run the heuristic returning either a Success wrapper around
        the result, or a Failure wrapper around the raised exceptions.
        """
        pass


class HFunction(Heuristic):
    """
    A wrapper around a unary function to make it act as a Heuristic.
    """
    def __init__(self, func):
        """Create a Heuristic from the given function."""
        self.func = func

    def __repr__(self):
        return 'HFunction(%r)' % self.func

    def __str__(self):
        return 'HFunction(%s)' % self.func

    def run(self, arg):
        """
        Run the wrapped function, returning a Result containing either
        the successful result value or the exceptions raised by the
        failure.
        """
        try:
            res = self.func(arg)
        except Exception as e:
            return Failure([e])
        else:
            return Success(res)


class HOrElses(Heuristic):
    """
    A composition of alternate Heuristics to calculate the same value
    different ways.

    If any succeeds, the composite succeeds.  If all fail, the
    composite fails.
    """

    def __init__(self, label, or_elses):
        """
        Create a Heuristic from an optional label and a list of the
        alternative Heuristics to try.
        """
        self.label = label
        self.or_elses = or_elses

    def __repr__(self):
        return 'HOrElses(%r, %r)' % (self.label, self.or_elses)

    def __str__(self):
        return 'HOrElses(%s, %s)' % (self.label, self.or_elses)

    def run(self, arg):
        """
        Run the alternative heuristics until one succeeds and return
        that value as a Success.  If none of them succeed, collect all
        their exceptions and return them as a Failure, labeling the
        group if the object has a label.
        """
        exceptions = []
        for h in self.or_elses:
            res = h.run(arg)
            if res.is_success():
                return res
            else:
                exceptions.extend(res.exceptions)

        # if you made it here, all of them failed
        if self.label:
            return Failure([(self.label, exceptions)])
        else:
            return Failure(exceptions)


class HAndThens(Heuristic):
    """
    A composition of Heuristics to be chained, each one's result
    becoming the input to the next.

    If all succeed, the composite succeeds.  If any fails, the
    composite fails.
    """

    def __init__(self, label, and_thens):
        """
        Create a Heuristic from an optional label and a list of
        Heuristics to be chained.
        """
        self.label = label
        self.and_thens = and_thens

    def run(self, arg):
        """
        Run the Heuristics chained, with the successful result value
        of each becoming the input argument to the next, and return
        that value as a Success.  If any of them fail, immediately
        abandon the rest of the calculation and return the raised
        exceptions as a Faillure, labeling the group if the object has
        a label.
        """
        if self.and_thens:
            for h in self.and_thens:
                res = h.run(arg)
                if res.is_success():
                    arg = res.value
                else:
                    if self.label:
                        return Failure([(self.label, res.exceptions)])
                    else:
                        return res
        else:
            res = Success(arg)

        # if you made it here, all of them succeeded
        return res


############################################################

def _toHeuristic(arg):
    """Wrap the argument if necessary and return an instance of Heuristic."""
    if isinstance(arg, Heuristic):
        return arg
    elif hasattr(arg, '__call__'):
        # it's a function, or something like one
        return HFunction(arg)
    else:
        raise ValueError('argument was neither a Heuristic nor callable: %r' %
                         arg)


def h_or_elses(*args):
    """
    Create an HOrElses using some syntactic sugar.  If the first
    argument is a string, it serves as a label; if not, the label is
    taken to be None.

    All other arguments become the components of the HOrElses.  The
    other arguments should all be Heuristics, but you may include
    unary functions, which will be wrapped to become HFunctions.
    """

    # Look for initial label
    if args and isinstance(args[0], str):
        label = args[0]
        args = args[1:]
    else:
        label = None

    return HOrElses(label, [_toHeuristic(arg) for arg in args])


def h_and_thens(*args):
    """
    Create an HAndThens using some syntactic sugar.  If the first
    argument is a string, it serves as a label; if not, the label is
    taken to be None.

    All other arguments become the components of the HAndThens.  The
    other arguments should all be Heuristics, but you may include
    unary functions, which will be wrapped to become HFunctions.
    """

    # Look for initial label
    if args and isinstance(args[0], str):
        label = args[0]
        args = args[1:]
    else:
        label = None

    return HAndThens(label, [_toHeuristic(arg) for arg in args])

############################################################


def _doubles(n): return n * 2


def _triples(n): return n * 3


def _succ(n): return n + 1


def _fails(n): raise Exception('I failed')


def _fails_with_message(msg):
    def _fails(arg):
        raise Exception(msg)
    return _fails


class TestHeuristic(unittest.TestCase):
    def test_is_success(self):
        self.assertTrue(Success(None).is_success())
        self.assertFalse(Failure([]).is_success())

    def test_is_failure(self):
        self.assertTrue(Failure([]).is_failure())
        self.assertFalse(Success(None).is_failure())

    def test_HFunction(self):
        res = HFunction(_doubles).run(3)
        self.assertTrue(res.is_success())
        self.assertEquals(6, res.value)

        res = HFunction(_fails).run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(1, len(res.exceptions))

    def test_h_or_elses(self):
        doubles = HFunction(_doubles)
        succ = HFunction(_succ)

        # Test contents
        h = h_or_elses(doubles, succ)
        self.assertTrue(isinstance(h, HOrElses))
        self.assertIsNone(h.label)
        self.assertEquals([doubles, succ], h.or_elses)

        h = h_or_elses('foo', succ, doubles)
        self.assertTrue(isinstance(h, HOrElses))
        self.assertEquals('foo', h.label)
        self.assertEquals([succ, doubles], h.or_elses)

        # Test function wrapping
        h = h_or_elses(_doubles, _succ)
        self.assertTrue(isinstance(h, HOrElses))
        self.assertIsNone(h.label)
        for o in h.or_elses:
            self.assertTrue(isinstance(o, HFunction))
        self.assertEquals([_doubles, _succ], [h.func for h in h.or_elses])

    def test_HOrElses(self):
        # I use h_or_elses() sugar here; the test above proves it to be
        # equivalent to explicit HOrElses construction.

        doubles = HFunction(_doubles)
        succ = HFunction(_succ)
        fails_foo = HFunction(_fails_with_message('foo'))
        fails_bar = HFunction(_fails_with_message('bar'))

        # Test that the first success returns, not trying further
        # alternatives
        res = h_or_elses(doubles, succ).run(3)
        self.assertTrue(res.is_success())
        self.assertEquals(6, res.value)

        res = h_or_elses(succ, doubles).run(3)
        self.assertTrue(res.is_success())
        self.assertEquals(4, res.value)

        res = h_or_elses(doubles, fails_foo).run(3)
        self.assertTrue(res.is_success())
        self.assertEquals(6, res.value)

        res = h_or_elses(fails_foo, doubles).run(3)
        self.assertTrue(res.is_success())
        self.assertEquals(6, res.value)

        # Test that multiple failures return all their exceptions.
        res = h_or_elses(fails_foo, fails_bar).run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(['foo', 'bar'],
                          [e.args[0] for e in res.exceptions])

        res = h_or_elses(fails_bar, fails_foo).run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(['bar', 'foo'],
                          [e.args[0] for e in res.exceptions])

        # Test that using a label wraps any exceptions.
        res = h_or_elses('wrapper', fails_bar, fails_foo).run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(1, len(res.exceptions))
        (lab, es) = res.exceptions[0]
        self.assertEquals('wrapper', lab)
        self.assertEquals(['bar', 'foo'], [e.args[0] for e in es])

        # Test empty case: no alternatives means immediate failure.
        res = h_or_elses().run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(0, len(res.exceptions))

    def test_h_and_thens(self):
        doubles = HFunction(_doubles)
        succ = HFunction(_succ)

        # Test contents
        h = h_and_thens(doubles, succ)
        self.assertTrue(isinstance(h, HAndThens))
        self.assertIsNone(h.label)
        self.assertEquals([doubles, succ], h.and_thens)

        h = h_and_thens('foo', succ, doubles)
        self.assertTrue(isinstance(h, HAndThens))
        self.assertEquals('foo', h.label)
        self.assertEquals([succ, doubles], h.and_thens)

        # Test function wrapping
        h = h_and_thens(_doubles, _succ)
        self.assertTrue(isinstance(h, HAndThens))
        self.assertIsNone(h.label)
        for a in h.and_thens:
            self.assertTrue(isinstance(a, HFunction))
        self.assertEquals([_doubles, _succ], [h.func for h in h.and_thens])

    def test_HAndThens(self):
        # I use h_and_thens() sugar here; the test above proves it to
        # be equivalent to explicit HAndThens construction.

        doubles = HFunction(_doubles)
        succ = HFunction(_succ)
        fails_foo = HFunction(_fails_with_message('foo'))
        fails_bar = HFunction(_fails_with_message('bar'))

        # Test that successful calculations run chained in order
        res = h_and_thens(doubles, succ).run(3)
        self.assertTrue(res.is_success())
        self.assertEquals(7, res.value)

        res = h_and_thens(succ, doubles).run(3)
        self.assertTrue(res.is_success())
        self.assertEquals(8, res.value)

        # Test that the first failure aborts the calculation
        res = h_and_thens(doubles, fails_foo).run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(1, len(res.exceptions))
        self.assertEquals('foo', res.exceptions[0].args[0])

        res = h_and_thens(fails_bar, doubles).run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(1, len(res.exceptions))
        self.assertEquals('bar', res.exceptions[0].args[0])

        res = h_and_thens(fails_foo, fails_bar).run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(1, len(res.exceptions))
        self.assertEquals('foo', res.exceptions[0].args[0])

        res = h_and_thens(fails_bar, fails_foo).run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(1, len(res.exceptions))
        self.assertEquals('bar', res.exceptions[0].args[0])

        # Test that using a label wraps any exceptions.
        res = h_and_thens('wrapper', fails_bar, fails_foo).run(3)
        self.assertFalse(res.is_success())
        self.assertEquals(1, len(res.exceptions))
        (lab, es) = res.exceptions[0]
        self.assertEquals('wrapper', lab)
        self.assertEquals(['bar'], [e.args[0] for e in es])

        # Test empty case: no actions means immediate success with no
        # change in argument.
        res = h_and_thens().run(3)
        self.assertTrue(res.is_success())
        self.assertEquals(3, res.value)


if __name__ == '__main__':
    unittest.main()
