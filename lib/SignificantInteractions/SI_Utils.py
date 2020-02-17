import pandas as pd
import os
import uuid
from ..installed_clients.WorkspaceClient import Workspace as Workspace
from ..installed_clients.DataFIleUtilClient import DataFileUtil


class SIintersect:
    """ Finds the intersection of significant interactions """

    def __init__(self, token, callback_url, scratch):
        self.token = token
        self.callback_url = callback_url
        self.scratch = scratch



