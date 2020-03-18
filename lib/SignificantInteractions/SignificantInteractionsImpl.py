# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid
import json
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
    def _df_to_list(self, df, fillna_val = 0, threshold=None):
        """
        _df_to_list: convert Dataframe to FloatMatrix2D matrix data
        """

        df.fillna(fillna_val, inplace=True)

        if threshold:
            drop_cols = list()
            for col in df.columns:
                if all(df[col] < threshold) and all(df[col] > -threshold):
                    drop_cols.append(col)
            df.drop(columns=drop_cols, inplace=True, errors='ignore')

            drop_idx = list()
            for idx in df.index:
                if all(df.loc[idx] < threshold) and all(df.loc[idx] > -threshold):
                    drop_idx.append(idx)
            df.drop(index=drop_idx, inplace=True, errors='ignore')

        matrix_data = {'row_ids': df.index.tolist(),
                       'col_ids': df.columns.tolist(),
                       'values': df.values.tolist()}

        return matrix_data

    def _save_corr_matrix(self, workspace_name, corr_matrix_name, corr_df, freq_df, sig_df, matrix_ref=None):
        """
        _save_corr_matrix: save
        KBaseExperiments.CorrelationMatrix
        object
        """

        if not isinstance(workspace_name, int):
            ws_name_id = self.dfu.ws_name_to_id(workspace_name)
        else:
            ws_name_id = workspace_name
        corr_data = {}
        corr_data.update({'coefficient_data': self._df_to_list(corr_df)})

        if matrix_ref:
            corr_data.update({'original_matrix_ref': matrix_ref})

        if sig_df is not None:
            corr_data.update({'significance_data': self._df_to_list(sig_df, 1)})

        if freq_df is not None:
            corr_data.update({'frequency_data': self._df_to_list(freq_df)})

        obj_type = 'KBaseExperiments.CorrelationMatrix'
        info = self.dfu.save_objects({
                "id": ws_name_id,
                "objects": [{
                    "type": obj_type,
                    "data": corr_data,
                    "name": corr_matrix_name
                    }]
                })[0]
        return "%s/%s/%s" % (info[6], info[0], info[4])
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.token = os.environ['KB_AUTH_TOKEN']
        self.wsURL = config['workspace-url']
        self.shared_folder = config['scratch']
        self.dfu = DataFileUtil(self.callback_url)
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

        logging.info('--->\nrunning metaMDS with input\n' +
                     'params:\n{}'.format(json.dumps(params, indent=1)))

        MatrixIds = params.get('MatrixIds')
        sig_cutoff = params.get('sig_cutoff')
        corr_cutoff = params.get('corr_cutoff')
        frequency = params.get('frequency')
        search_for_type = params.get('search_for_type')
        if sig_cutoff is None and corr_cutoff is None:
            raise ValueError("ERROR: Both sig_cutoff and corr_cutoff are null. At least one is needed")
        if frequency is None:
            frequency = 0
        corr_matrix_name = params.get('corr_matrix_name')

        si = SI(token=self.token, callback_url=self.callback_url, scratch=self.shared_folder)
        si_dict = si.run(MatrixIds=MatrixIds, sig_cutoff=sig_cutoff, corr_cutoff=corr_cutoff, frequency=frequency,
                         search_for_type=search_for_type)

        corr_matrix_obj_ref = self._save_corr_matrix(workspace_name=params['workspace_name'],
                                                     corr_matrix_name=corr_matrix_name,
                                                     corr_df=si_dict['corr_df'], sig_df=si_dict['sig_df'],
                                                     freq_df=si_dict['freq_df'], matrix_ref=MatrixIds)

        report_client = KBaseReport(self.callback_url, token=self.token)
        report_name = "Significant_Interaction_Intersect_" + str(uuid.uuid4())
        report_info = report_client.create_extended_report({
            'objects_created': [{'ref': corr_matrix_obj_ref,
                                 'description': 'Correlation Matrix'}],
            'direct_html_link_index': 0,
            'html_links': si_dict['html_paths'],
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
