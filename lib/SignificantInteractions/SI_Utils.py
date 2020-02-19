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
        self.a_dict = {}
        self.html_paths = []
        self.corr_df = None
        self.sig_df = None
        self.freq_df = None
        self.dfu = DataFileUtil(self.callback_url)

    # Returns Significance Matrix pd.DataFrame(): get_pd_matrix(MatrixId)
    def get_pd_matrix(self, MatrixId):

        obj = self.dfu.get_objects({'object_refs': [MatrixId]})

        co = obj['data'][0]['data']['coefficient_data']
        co_rows = co['row_ids']
        co_cols = co['col_ids']
        co_vals = co['values']

        co_mat = pd.DataFrame(co_vals, index=co_rows, columns=co_cols)

        sig = obj['data'][0]['data']['significance_data']
        sig_rows = sig['row_ids']
        sig_cols = sig['col_ids']
        sig_vals = sig['values']

        sig_mat = pd.DataFrame(sig_vals, index=sig_rows, columns=sig_cols)

        return [sig_mat, co_mat]

    #
    def push_to_dict(self, sig_matrix, co_matrix, sig_cutoff, corr_cutoff):
        otu_1s = sig_matrix.index
        otu_2s = sig_matrix.columns
        for i in range(len(sig_matrix.index)):
            for j in range(i + 1, len(sig_matrix.index)):
                key = otu_1s[i] + '<->' + otu_2s[j]
                sig_val = sig_matrix.iloc[i][j]
                co_val = co_matrix[otu_1s[i]][otu_2s[j]]
                if sig_val <= sig_cutoff and co_val >= corr_cutoff:
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

        # lists
        row_col_list = []
        val_list = []
        # Make dict to make html file
        html_dict = {}
        for key, val in self.a_dict.items():
            # test criteria for html_dict and DataFrame
            if val[2] >= frequency:
                html_dict.update({key: [val[0] / quantity, val[1] / quantity, val[2]]})
                OTUs = key.split('<->')
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
            mats = self.get_pd_matrix(MatrixId=Id)
            self.push_to_dict(sig_matrix=mats[0], co_matrix=mats[1], sig_cutoff=sig_cutoff, corr_cutoff=corr_cutoff)
        self.to_html(frequency=frequency, quantity=len(MatrixIds))
        return {
            'html_paths': self.html_paths,
            'corr_df': self.corr_df,
            'sig_df': self.sig_df,
            'freq_df': self.freq_df
        }


"""class SIintersect:
    ''' Finds the intersection of significant interactions '''

    def __init__(self, token, callback_url, scratch):
        self.token = token
        self.callback_url = callback_url
        self.scratch = scratch
        self.each_intersect_dict = {}

    # Returns Significance Matrix pd.DataFrame(): get_pd_matrix(MatrixId)
    def get_pd_matrix(self, MatrixId):

        dfu = DataFileUtil(self.callback_url)
        obj = dfu.get_objects({'object_refs': [MatrixId]})

        sig = obj['data'][0]['data']['significance_data']
        sig_rows = sig['row_ids']
        sig_cols = sig['col_ids']
        sig_vals = sig['values']

        mat = pd.DataFrame(index=sig_rows, columns=sig_cols)
        for i in range(len(mat.index)):
            mat.iloc[i] = sig_vals[i]
        return mat

    # Returns dictionary of significant interactions based on cutoff: get_significant_dict(matrix, cutoff=0.7)
    def get_significant_dict(self, matrix, cutoff=0.7):
        otu_1s = matrix.index
        otu_2s = matrix.columns
        the_dict = {}
        for i in range(len(matrix.index)):
            for j in range(i + 1, len(matrix.index)):
                key = otu_1s[i] + '<->' + otu_2s[j]
                val = matrix.iloc[i][j]
                if val >= cutoff:
                    try:
                        the_dict[key] += val
                    except KeyError:
                        the_dict.update({key: val})
        return the_dict

    # Returns dictionary of elements that intersect: get_intersection_dict(dict1={}, dict2={})
    def get_intersection_dict(self, dict1, dict2):
        intersection_dict = {}

        for key, val in dict2.items():
            if key in dict1.keys():
                intersection_dict.update({key: val + dict1[key]})
        return intersection_dict

    # Make dictionary of dictionaries of each intersection comparing neighboring matrix--(1,2), (2,3), ..., (n,1).
    def get_dict_of_each(self, MatrixIds, cutoff):
        key1 = MatrixIds[-1]
        m1 = self.get_pd_matrix(MatrixId=MatrixIds[-1])
        dict1 = self.get_significant_dict(matrix=m1, cutoff=cutoff)
        for i in range(len(MatrixIds)):
            key2 = MatrixIds[i]
            m2 = self.get_pd_matrix(MatrixId=MatrixIds[i])
            dict2 = self.get_significant_dict(matrix=m2, cutoff=cutoff)
            the_dict = self.get_intersection_dict(dict1=dict1, dict2=dict2)
            self.each_intersect_dict.update({key1 + '<->' + key2: the_dict})
            key1 = key2
            dict1 = dict2
        return self.each_intersect_dict

    # run
    def run(self, MatrixIds, cutoff):
        m = self.get_pd_matrix(MatrixIds[0])
        main_dict = self.get_significant_dict(m, cutoff)
        for Id in MatrixIds[1:]:
            m = self.get_pd_matrix(Id)
            dict2 = self.get_significant_dict(m, cutoff)
            main_dict = self.get_intersection_dict(dict1=main_dict, dict2=dict2)
        for key, val in main_dict.items():
            main_dict[key] = val / len(MatrixIds)
        return main_dict"""