import datetime
import pickle

class Report:
    def __init__(self, name="noname"):
        self.name = str(datetime.datetime.now()) + "_" + name
        self.dictionary = {
            "name"  : self.name
        }

    def add_information(self, key, data):
        self.dictionary[key] = data

    def save(self):
        with open(self.name + ".pickle", mode="wb") as f:
            pickle.dump(self.dictionary, f)