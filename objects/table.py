class Job:
    def __init__(self, conditions):
        self.result     = None
        self.end_flag   = False
        self.__dict__.update(conditions)

class JobTable:
    def __init__(self):
        self.reset()

    def submit(self, job):
        self.table.append(job)

    def reset(self):
        self.table  = []
