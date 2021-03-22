from pymechlib.cmechlib import _FeedForward, _Pid
from pymechlib import rnd
import numpy as np


class FeedForward(_FeedForward):

    def __init__(self, period, amplitudes, phases):

        self.amplitudes = np.array(amplitudes.copy(), dtype=np.float64)
        self.phases = np.array(phases.copy(), dtype=np.float64)
        self.period = period
        self.order = len(amplitudes)

        _FeedForward.__init__(self, 2*np.pi/(3*period), amplitudes, phases)

    def mutate(self):

        new_amplitudes = self.amplitudes
        new_phases = self.phases

        # Modify the parameters of the controller randomly
        pick = rnd.choice([0, 1, 2])

        if pick == 0:
            new_phases += rnd.uniform(-0.1, 0.1)
        if pick == 1:
            new_amplitudes += rnd.uniform(-0.1, 0.1)
        else:
            new_amplitudes = new_amplitudes * rnd.uniform(0.9, 1.1)

        return FeedForward(self.period, new_amplitudes, new_phases)

    def copy(self):

        return FeedForward(self.period, self.amplitudes, self.phases)


class Pid(_Pid):

    def __init__(self, time, parameters):

        self.parameters = parameters.copy()
        self.time = time

        _Pid.__init__(self, time, parameters)

    def mutate(self):

        new_parameters = self.parameters.copy()

        pick = rnd.choice(3)

        if pick == 0:
            new_parameters[0] += rnd.uniform(-0.2, 0.2)
        elif pick == 1:
            new_parameters[1] += rnd.uniform(-0.1, 0.1)
        else:
            new_parameters[2] += rnd.uniform(-0.05, 0.05)

        return Pid(self.time, new_parameters)

    def copy(self):

        return Pid(self.time, self.parameters)

