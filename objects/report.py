import datetime
import pickle
import os

REPORT_DIRECTORY = 'C:\\Users\\user\\Documents\\scalabrad\\heya\\reports\\'

class Report:
    def __init__(self, name="noname"):
        now = datetime.datetime.now()
        self.name = now.strftime("%Y_%m%d_%H%M%S") + "_" + name
        self.dictionary = {
            "name"  : self.name
        }

    def add_information(self, key, data):
        self.dictionary[key] = data

    def save(self):
        with open(REPORT_DIRECTORY + self.name + ".pickle", mode="wb") as f:
            pickle.dump(self.dictionary, f)