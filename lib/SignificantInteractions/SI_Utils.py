import pandas as pd
import os
import uuid
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.DataFIleUtilClient import DataFileUtil


# Class2
class SI:

    def __init__(self, token, callback_url, scratch):
        self.token = token
        self.callback_url = callback_url
        self.scratch = scratch
        self.a_dict = {}  # values are list in form [sig_val, corr_val, frequency]
        self.html_paths = []
        self.corr_df = None
        self.sig_df = None
        self.freq_df = None
        self.dfu = DataFileUtil(self.callback_url)

    # Returns Correlation and Significance Matrix pd.DataFrame()
    def get_pd_matrix(self, MatrixId, corr_cutoff, sig_cutoff):
        # Initialize dictionary to be returned
        returning_dict = {
            'corr_mat': None,
            'sig_mat': None
        }
        obj = self.dfu.get_objects({'object_refs': [MatrixId]})

        # If 'coefficient_data' exist
        if obj['data'][0]['data']['coefficient_data']:
            co = obj['data'][0]['data']['coefficient_data']
            co_rows = co['row_ids']
            co_cols = co['col_ids']
            co_vals = co['values']
            co_mat = pd.DataFrame(co_vals, index=co_rows, columns=co_cols)
            returning_dict['corr_mat'] = co_mat

        # If 'significant_data' exist
        if obj['data'][0]['data']['significance_data']:
            sig = obj['data'][0]['data']['significance_data']
            sig_rows = sig['row_ids']
            sig_cols = sig['col_ids']
            sig_vals = sig['values']
            sig_mat = pd.DataFrame(sig_vals, index=sig_rows, columns=sig_cols)
            returning_dict['sig_mat'] = sig_mat

        return returning_dict

    # Push matrix keys(index<->column) and values into dictionary
    def push_to_dict(self, matrix_dict, sig_cutoff, corr_cutoff):
        if sig_cutoff is not None and corr_cutoff is not None:
            otu_1s = matrix_dict['sig_mat'].index
            otu_2s = matrix_dict['sig_mat'].columns
            # Go through half of matrix since keys/values repeat, keys just in reverse
            for i in range(len(matrix_dict['sig_mat'].index)):
                for j in range(i + 1, len(matrix_dict['sig_mat'].index)):
                    # Get otu's for key and then sort to cover both possibilities
                    # so won't be seperate when pushing into dict
                    sorted_otus = [otu_1s[i], otu_2s[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    sig_val = matrix_dict['sig_mat'][otu_1s[i]][otu_2s[j]]
                    co_val = matrix_dict['corr_mat'][otu_1s[i]][otu_2s[j]]
                    # Increment frequency if
                    if sig_val <= sig_cutoff and co_val >= corr_cutoff:
                        try:
                            self.a_dict[key][0] += sig_val
                            self.a_dict[key][1] += co_val
                            self.a_dict[key][2] += 1
                        except KeyError:
                            self.a_dict.update({key: [sig_val, co_val, 1]})
                    # else just do this
                    else:
                        try:
                            self.a_dict[key][0] += sig_val
                            self.a_dict[key][1] += co_val
                        except KeyError:
                            self.a_dict.update({key: [sig_val, co_val, 0]})

        # If no sig_cutoff is specified
        elif sig_cutoff is not None:
            otu_1s = matrix_dict['sig_mat'].index
            otu_2s = matrix_dict['sig_mat'].columns
            for i in range(len(matrix_dict['sig_mat'].index)):
                for j in range(i + 1, len(matrix_dict['sig_mat'].index)):
                    sorted_otus = [otu_1s[i], otu_2s[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    if matrix_dict['corr_mat'] is not None:
                        co_val = matrix_dict['corr_mat'][otu_1s[i]][otu_2s[j]]
                    else:
                        co_val = 0
                    sig_val = matrix_dict['sig_mat'][otu_1s[i]][otu_2s[j]]
                    if sig_val <= sig_cutoff:
                        try:
                            self.a_dict[key][0] += sig_val
                            self.a_dict[key][1] += co_val
                            self.a_dict[key][2] += 1
                        except KeyError:
                            self.a_dict.update({key: [sig_val, co_val, 1]})
                    else:
                        try:
                            self.a_dict[key][0] += sig_val
                            self.a_dict[key][1] += co_val
                        except KeyError:
                            self.a_dict.update({key: [sig_val, co_val, 0]})

        # if no corr_cutoff is specified
        elif corr_cutoff is not None:
            otu_1s = matrix_dict['corr_mat'].index
            otu_2s = matrix_dict['corr_mat'].columns
            for i in range(len(matrix_dict['corr_mat'].index)):
                for j in range(i + 1, len(matrix_dict['corr_mat'].index)):
                    sorted_otus = [otu_1s[i], otu_2s[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    if matrix_dict['sig_mat'] is not None:
                        sig_val = matrix_dict['sig_mat'][otu_1s[i]][otu_2s[j]]
                    else:
                        sig_val = 0
                    co_val = matrix_dict['corr_mat'][otu_1s[i]][otu_2s[j]]
                    if co_val >= corr_cutoff:
                        try:
                            self.a_dict[key][0] += sig_val
                            self.a_dict[key][1] += co_val
                            self.a_dict[key][2] += 1
                        except KeyError:
                            self.a_dict.update({key: [sig_val, co_val, 1]})
                    else:
                        try:
                            self.a_dict[key][0] += sig_val
                            self.a_dict[key][1] += co_val
                        except KeyError:
                            self.a_dict.update({key: [sig_val, co_val, 0]})

    def to_html(self, frequency, quantity):
        # set up directory in scratch
        output_dir = os.path.join(self.scratch, str(uuid.uuid4()))
        os.mkdir(output_dir)
        # set up directory for html folder
        html_folder = os.path.join(output_dir, 'html')
        os.mkdir(html_folder)

        # list for index and columns
        row_col_list = []
        # Make dict to make html file
        html_dict = {}
        for key, val in self.a_dict.items():
            # test criteria for html_dict and DataFrame
            if val[2] >= frequency:
                html_dict.update({key: [val[0] / quantity, val[1] / quantity, val[2]]})
                OTUs = key.split('<->')
                # add otu's to row_col_list
                if OTUs[0] not in row_col_list:
                    row_col_list.append(OTUs[0])
                if OTUs[1] not in row_col_list:
                    row_col_list.append(OTUs[1])
        # sort row_col_list
        row_col_list.sort()
        # pandas DataFrame
        self.corr_df = pd.DataFrame(index=row_col_list, columns=row_col_list)
        self.sig_df = pd.DataFrame(index=row_col_list, columns=row_col_list)
        self.freq_df = pd.DataFrame(index=row_col_list, columns=row_col_list)
        # Make html_str out of html_dict
        html_str = "<html>" \
                   "<body>" \
                   '<table border="2">' \
                   "<tr>" \
                   "<td>OTUs: </td><td>Average Significance: </td> <td>Average Correlation: </td> <td>Frequency:</td>" \
                   "</tr>"
        for key, val in html_dict.items():
            # Push values into df matrices
            OTUs = key.split('<->')
            self.corr_df[OTUs[0]][OTUs[1]] = val[1]
            self.corr_df[OTUs[1]][OTUs[0]] = val[1]
            self.sig_df[OTUs[0]][OTUs[1]] = val[0]
            self.sig_df[OTUs[1]][OTUs[0]] = val[0]
            self.freq_df[OTUs[0]][OTUs[1]] = val[2]
            self.freq_df[OTUs[1]][OTUs[0]] = val[2]
            # html part
            html_str += "<tr>" \
                        "<td>" + key + ":</td><td>" + str(round(val[0], 5)) + "</td><td>" + str(round(val[1], 5)) \
                        + "</td><td>" + str(val[2]) + " / " + str(quantity) + "</td>" \
                            "</tr>"
        html_str += "</table>" \
                    "</body>" \
                    "</html>"

        with open(os.path.join(html_folder, "index.html"), 'w') as index_file:
            index_file.write(html_str)

        # have needed files saved to folder before shock
        shock = self.dfu.file_to_shock({'file_path': html_folder,
                                        'make_handle': 0,
                                        'pack': 'zip'})
        # list that goes to 'html_links'
        self.html_paths.append({'shock_id': shock['shock_id'],
                                'name': 'index.html',
                                'label': 'html files',
                                'description': "desc"})

    def run(self, MatrixIds, sig_cutoff, corr_cutoff, frequency):
        for Id in MatrixIds:
            mats = self.get_pd_matrix(MatrixId=Id, corr_cutoff=corr_cutoff, sig_cutoff=sig_cutoff)
            self.push_to_dict(matrix_dict=mats, sig_cutoff=sig_cutoff, corr_cutoff=corr_cutoff)
        self.to_html(frequency=frequency, quantity=len(MatrixIds))
        return {
            'html_paths': self.html_paths,
            'corr_df': self.corr_df,
            'sig_df': self.sig_df,
            'freq_df': self.freq_df
        }
