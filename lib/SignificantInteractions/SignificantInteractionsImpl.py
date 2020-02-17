# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
from SignificantInteractions.SI_Utils import SI
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFIleUtilClient import DataFileUtil
#END_HEADER


class SignificantInteractions:
    '''
    Module Name:
    SignificantInteractions

    Module Description:
    A KBase module: SignificantInteractions
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.token = os.environ['KB_AUTH_TOKEN']
        self.wsURL = config['workspace-url']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_SignificantInteractions(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_SignificantInteractions

        MatrixIds = params.get('MatrixIds')
        cutoff = params.get('cutoff')
        frequency = params.get('frequency')

        si = SI(token=self.token, callback_url=self.callback_url, scratch=self.shared_folder)
        html_paths = si.run(MatrixIds=MatrixIds, cutoff=cutoff, frequency=frequency)

        report_client = KBaseReport(self.callback_url, token=self.token)
        report_name = "Significant_Interaction_Intersect_" + str(uuid.uuid4())
        report_info = report_client.create_extended_report({
            'direct_html_link_index': 0,
            'html_links': html_paths,
            'report_object_name': report_name,
            'workspace_name': params['workspace_name']
        })
        output = {
            'report_ref': report_info['ref'],
            'report_name': report_info['name'],
        }
        #END run_SignificantInteractions

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_SignificantInteractions return value ' +
                             'output is not type dict as required.')
        # return the results
        return [report_info]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
