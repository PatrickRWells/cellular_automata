import random

import matplotlib.pyplot as plt
from numpy import array_equal


def random_string(length):

    """
    Returns a random string of a given length, where each character
    is a 0, 1, or 2.

    Parameters
    ----------
    length: int
        Posivite integer that specifies the desired length of the string.

    Returns
    -------
    out: list
        The random  string given as a list, with int elements.
    """

    if not isinstance(length, int) or length < 0:
        raise ValueError("input length must be a positive ingeter")
    return [random.randint(0,2) for _ in range(length)]


class TCA:
    def __init__(self, rule_number = 0):

        '''
        Initializes a tri-state automata simulator
        Parameters
        ----------
        Rule number: int
            Integer value between 0 and 19682, inclusive

        Attributes
        ----------
        rule_number: int
            Integer value between 0 and 19682, inclusive.
            Specifies the rules under which the automata evolves
        rule_number_trinary: string
            rule_number represented in base 3
        lookup:
            A lookup table for evolving the cellular automata
            as specified by the rule number
        '''

        self.rule_number = rule_number
        self.rule_number_trinary = self.to_trinary(rule_number)
        self.lookup = self.lookup_table()

    def update_rule(self,rule_number):

        """
        Updates a simulator to use a new rule

        Parameters
        ----------
        Rule number: int
            Integer value between 0 and 19682, inclusinve

        Returns
        -------
        Nothing. State of the object is updated
        """

        self.rule_number = rule_number
        self.rule_number_trinary = self.to_trinary(rule_number)
        self.lookup = self.lookup_table()


    def to_trinary(self, rule_number):

        """
        Returns the value of the input in base-3

        Parameters
        ----------
        rule_number: int
            Integer value between 0 and 19682, inclusive.

        Returns
        -------
        trinary: string
            Input number in base 3, as a string
        """

        if (not isinstance(rule_number, int)
           or rule_number < 0
           or rule_number > 19682):
            raise ValueError("rule_number must be an int " \
                             "between 0  and 19682, inclusive")

        in_trinary = ''
        tmp = rule_number
        for i in range(8, -1, -1):
            factor = int(pow(3,i))
            val = tmp // factor
            in_trinary += str(val)
            tmp -= factor*val
        return(in_trinary)

    def lookup_table(self):

        """
        Returns a dictionary which maps ECA neighborhoods to output values.
        Uses Wolfram rule number convention.

        Parameters
        ----------
        None

        Returns
        -------
        lookup_table: dict
            Lookup table dictionary that maps neighborhood tuples to
            their output according to the ECA local evolution rule
            (i.e. the lookup table), as specified by the rule number.
        """

        neighborhoods = [(0,0), (0,1), (0,2), (1,0), (1,1),
                         (1,2), (2,0), (2,1), (2,2)]

        return dict(zip(neighborhoods,
                    map(int,reversed(self.rule_number_trinary))))
        # use map so that outputs are ints, not strings


    def spacetime_field(self, initial_condition, time_steps):

        """

        Returns a spacetime field array using the given rule number on the
        given initial condition for the given number of time steps.

        Parameters
        ----------
        initial_condition: list
            Trinary string used as the initial condition for the ECA.
            Elements of the list should be ints.
        time_steps: int
            Positive integer specifying
            the number of time steps for evolving the ECA.

        Returns:
            Spacetime field consisting of initial conditions evolved
            for the given number of time steps

        """


        if time_steps < 0:
            raise ValueError("time_steps must be a non-negative integer")
        try:
            time_steps = int(time_steps)
        except ValueError:
            raise ValueError("time_steps must be a non-negative integer")

        for i in initial_condition:
            if i not in [0,1,2]:
                raise ValueError("initial condition must be a " \
                                "list of 0s, 1s, and 2s")

        # initialize spacetime field and current configuration
        spacetime_field = [initial_condition]
        current_configuration = initial_condition.copy()
        length = len(current_configuration)

        # apply the lookup table to evolve the CA
        # for the given number of time steps

        for t in range(time_steps):
            new_configuration = []

            for i in range(length):
                neighborhood = (current_configuration[(i-1)],
                                current_configuration[i])
                new_configuration.append(self.lookup[neighborhood])

            current_configuration = new_configuration
            spacetime_field.append(new_configuration)

        return spacetime_field

    def spacetime_diagram(self, spacetime_field, size=12, colors=plt.cm.Greys):

        """
        Produces a simple spacetime diagram image using matplotlib
        imshow with 'nearest' interpolation.

       Parameters
        ---------
        spacetime_field: array-like (2D)
            1+1 dimensional spacetime field, given as a 2D array
            or list of lists. Time should be dimension 0 so that
            spacetime_field[t] is the spatial configuration at time t.

        size: int, optional (default=12)
            Sets the size of the figure: figsize=(size,size)
        colors: matplotlib colormap, optional (default=plt.cm.Greys)
            See https://matplotlib.org/tutorials/colors/colormaps.html
            for colormap choices. A colormap 'cmap' is called as
            colors=plt.cm.cmap

        """

        plt.figure(figsize=(size,size))
        plt.imshow(spacetime_field, cmap=colors, interpolation='nearest', origin='lower')
        plt.show()

    def simulate(self, initial_conditions, time_steps, figsize=12):

        """
        Run the cellular automata for a given number of time steps starting
        from some initial conditions

        Parameters
        ----------
        initial_conditions: list
            Trinary string used as the initial condition for the ECA.
            Elements of the list should be ints.
        time_steps: int
            Positive integer indicating number of time steps

        figsize: int
            size of the resultant figure

        """

        field = self.spacetime_field(initial_conditions, time_steps)
        self.spacetime_diagram(field, figsize)

def test_rule0():
    tca = TCA()
    expected_out = [0 for i in range(10)]
    initial = random_string(10)
    field = tca.spacetime_field(initial, 2)
    
    for t, config in enumerate(field[1:]):
        assert array_equal(config, expected_out), \
        "Rule 0 test failed. Configuration after {} steps incorrect".format(t)

    print("Rule 0 test passed")

def test_rule9841():
    tca = TCA(9841)
    expected_out = [1 for i in range(10)]
    initial = random_string(10)
    field = tca.spacetime_field(initial, 2)

    for t, config in enumerate(field[1:]):
        assert array_equal(config, expected_out), \
        "Rule 9841 test failed. Configuration after {} steps incorrect".format(t)

    print("Rule 9841 test passed")

def test_rule19682():
    tca = TCA(19682)
    expected_out = [2 for i in range(10)]
    initial = random_string(10)
    field = tca.spacetime_field(initial, 2)

    for t, config in enumerate(field[1:]):
        assert array_equal(config, expected_out), \
        "Rule 19682 test failed. Configuration after {} steps incorrect".format(t)

    print("Rule 19682 test passed")




if __name__ == "__main__":
    test_rule0()
    test_rule9841()
    test_rule19682()