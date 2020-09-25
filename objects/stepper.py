DEBUG_MODE = False

def set_debug_mode(mode):
    global DEBUG_MODE
    DEBUG_MODE = mode

class Stepper:
    def __init__(self,execute,n_param):
        self.execute    = execute
        self.n_param    = n_paraｍ
        self.reset()

    def step(self,phi,evaluate=False,stop_signal=False):
        report          = self.execute(phi)
        score           = report["score"]

        if evaluate:
            self.score.append(score)
            self.phi.append(phi)
            self.iteration.append(self.step_count)
            self.register.append(report["register"])
            self.report(report)
        else:
            self.step_count += report["step"]

        return score

    def callback(self,phi):
        self.step(phi,evaluate=True)

    def reset(self):
        self.step_count = 0
        self.score      = []
        self.phi        = []
        self.iteration  = []
        self.register   = []

    def report(self,report):
        print('\033[34m')
        print("Iteration : [{0}]".format(len(self.iteration)))
        if DEBUG_MODE:
            print(" "*5 + str(report))
        else:
            print(" "*5 + str(report["score"]))
        print('\033[0m')