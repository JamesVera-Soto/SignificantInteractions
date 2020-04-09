import pandas as pd
import numpy as np
import os
import uuid
import logging
from installed_clients.WorkspaceClient import Workspace as Workspace
from installed_clients.DataFIleUtilClient import DataFileUtil


# Class2
class SI:

    def __init__(self, token, callback_url, scratch):
        self.token = token
        self.callback_url = callback_url
        self.scratch = scratch
        self.a_dict = {}  # values are list in form [corr_val, sig_val, frequency]
        self.html_paths = []
        self.corr_df = None
        self.sig_df = None
        self.freq_df = None
        self.dfu = DataFileUtil(self.callback_url)
        self.sig_cutoff = None
        self.corr_cutoff = None
        # Object info
        self.obj = None
        self.sig_rows = None
        self.sig_cols = None
        self.sig_vals = None
        self.corr_rows = None
        self.corr_cols = None
        self.corr_vals = None
        self.is_unique_search = False

    # Returns Correlation and Significance Matrix pd.DataFrame()
    def _get_matrix_obj(self, MatrixId, corr_cutoff, sig_cutoff):

        logging.info('getting matrix: {}'.format(MatrixId))

        self.obj = self.dfu.get_objects({'object_refs': [MatrixId]})

        # If 'coefficient_data' exist
        try:
            corr = self.obj['data'][0]['data']['coefficient_data']
            self.corr_rows = corr['row_ids']
            self.corr_cols = corr['col_ids']
            self.corr_vals = corr['values']
            self.corr_cutoff = corr_cutoff
        except KeyError:
            pass

        # If 'significant_data' exist
        try:
            sig = self.obj['data'][0]['data']['significance_data']
            self.sig_rows = sig['row_ids']
            self.sig_cols = sig['col_ids']
            self.sig_vals = sig['values']
            self.sig_cutoff = sig_cutoff
        except KeyError:
            pass

    # Push matrix keys(index<->column) and values into dictionary
    def _push_to_dict(self, sig_cutoff, corr_cutoff):
        """
        Puts all interactions into a dictionary
        """

        logging.info('push_to_dict with corr_cutoff: {} , and sig_cutoff: {}'.format(corr_cutoff, sig_cutoff))

        length = len(self.corr_rows)
        # If both sig_data and corr_data are used
        if sig_cutoff is not None and corr_cutoff is not None:
            # Go through half of matrix since keys/values repeat, keys just in reverse
            for i in range(length):
                if i % 1000 == 0:
                    logging.info('updating... on row# ' + str(i))
                for j in range(i + 1, length):
                    # Get otu's for key and then sort to cover both possibilities
                    # so won't be seperate when pushing into dict
                    sorted_otus = [self.corr_rows[i], self.corr_cols[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    sig_val = self.sig_vals[i][j]
                    corr_val = self.corr_vals[i][j]
                    # Increment frequency if
                    if sig_val <= sig_cutoff and corr_val >= corr_cutoff:
                        try:
                            self.a_dict[key][0] += corr_val
                            self.a_dict[key][1] += sig_val
                            self.a_dict[key][2] += 1
                        except KeyError:
                            self.a_dict.update({key: [corr_val, sig_val, 1]})
                    # else just do this
                    else:
                        try:
                            self.a_dict[key][0] += corr_val
                            self.a_dict[key][1] += sig_val
                        except KeyError:
                            self.a_dict.update({key: [corr_val, sig_val, 0]})

        # If no corr_cutoff is specified but sig_cutoff is; or corr_matrix doesn't exit
        elif sig_cutoff is not None:
            for i in range(length):
                for j in range(i + 1, length):
                    if i % 1000 == 0:
                        logging.info('updating... on row# ' + str(i))
                    sorted_otus = [self.sig_rows[i], self.sig_cols[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    if self.corr_cutoff is not None:
                        corr_val = self.corr_vals[i][j]
                    else:
                        corr_val = 0
                    sig_val = self.sig_vals[i][j]
                    if sig_val <= sig_cutoff:
                        try:
                            self.a_dict[key][0] += corr_val
                            self.a_dict[key][1] += sig_val
                            self.a_dict[key][2] += 1
                        except KeyError:
                            self.a_dict.update({key: [corr_val, sig_val, 1]})
                    else:
                        try:
                            self.a_dict[key][0] += corr_val
                            self.a_dict[key][1] += sig_val
                        except KeyError:
                            self.a_dict.update({key: [corr_val, sig_val, 0]})

        # if no sig_cutoff is specified but corr_cutoff is; or sig_matrix doesn't exist
        elif corr_cutoff is not None:
            for i in range(length):
                if i % 1000 == 0:
                    logging.info('updating... on row# ' + str(i))
                for j in range(i + 1, length):
                    sorted_otus = [self.corr_rows[i], self.corr_cols[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    if self.sig_cutoff is not None:
                        sig_val = self.sig_vals[i][j]
                    else:
                        sig_val = 0
                    corr_val = self.corr_vals[i][j]
                    if corr_val >= corr_cutoff:
                        try:
                            self.a_dict[key][0] += corr_val
                            self.a_dict[key][1] += sig_val
                            self.a_dict[key][2] += 1
                        except KeyError:
                            self.a_dict.update({key: [corr_val, sig_val, 1]})
                    else:
                        try:
                            self.a_dict[key][0] += corr_val
                            self.a_dict[key][1] += sig_val
                        except KeyError:
                            self.a_dict.update({key: [corr_val, sig_val, 0]})
        else:
            raise ValueError('ERROR: no comparing can be performed. Perhaps no corr_cutoff was specified and '
                             'significance data does not exist.')

    def _push_to_unique_dict(self, sig_cutoff, corr_cutoff):
        """
        Puts the interactions that meet the criteria of the specified unique matrix into a dictionary
        """

        logging.info('_push_to_unique_dict with corr_cutoff: {} , and sig_cutoff: {}'.format(corr_cutoff, sig_cutoff))

        self.is_unique_search = True
        length = len(self.corr_rows)
        # If both sig_data and corr_data are used
        if sig_cutoff is not None and corr_cutoff is not None:
            # Go through half of matrix since keys/values repeat, keys just in reverse
            for i in range(length):
                if i % 1000 == 0:
                    logging.info('updating... on row# ' + str(i))
                for j in range(i + 1, length):
                    # Get otu's for key and then sort to cover both possibilities
                    # so won't be seperate when pushing into dict
                    sorted_otus = [self.corr_rows[i], self.corr_cols[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    sig_val = self.sig_vals[i][j]
                    corr_val = self.corr_vals[i][j]
                    # Increment frequency if
                    if sig_val <= sig_cutoff and corr_val >= corr_cutoff:
                        try:
                            self.a_dict[key][0] += corr_val
                            self.a_dict[key][1] += sig_val
                            self.a_dict[key][2] += 1
                        except KeyError:
                            self.a_dict.update({key: [corr_val, sig_val, 1]})

        # If no corr_cutoff is specified but sig_cutoff is; or corr_matrix doesn't exit
        elif sig_cutoff is not None:
            for i in range(length):
                for j in range(i + 1, length):
                    if i % 1000 == 0:
                        logging.info('updating... on row# ' + str(i))
                    sorted_otus = [self.sig_rows[i], self.sig_cols[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    if self.corr_cutoff is not None:
                        corr_val = self.corr_vals[i][j]
                    else:
                        corr_val = 0
                    sig_val = self.sig_vals[i][j]
                    if sig_val <= sig_cutoff:
                        try:
                            self.a_dict[key][0] += corr_val
                            self.a_dict[key][1] += sig_val
                            self.a_dict[key][2] += 1
                        except KeyError:
                            self.a_dict.update({key: [corr_val, sig_val, 1]})

        # if no sig_cutoff is specified but corr_cutoff is; or sig_matrix doesn't exist
        elif corr_cutoff is not None:
            for i in range(length):
                if i % 1000 == 0:
                    logging.info('updating... on row# ' + str(i))
                for j in range(i + 1, length):
                    sorted_otus = [self.corr_rows[i], self.corr_cols[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    if self.sig_cutoff is not None:
                        sig_val = self.sig_vals[i][j]
                    else:
                        sig_val = 0
                    corr_val = self.corr_vals[i][j]
                    if corr_val >= corr_cutoff:
                        try:
                            self.a_dict[key][0] += corr_val
                            self.a_dict[key][1] += sig_val
                            self.a_dict[key][2] += 1
                        except KeyError:
                            self.a_dict.update({key: [corr_val, sig_val, 1]})
        else:
            raise ValueError('ERROR: no comparing can be performed. Perhaps no corr_cutoff was specified and '
                             'significance data does not exist.')

    def _remove_from_unique_dict(self, sig_cutoff, corr_cutoff):
        """
        Removes interactions from dictionary if they are found meeting the criteria in another matrix
        """

        logging.info('remove_from_unique_dict with corr_cutoff: {} , and sig_cutoff: {}'.format(corr_cutoff,
                                                                                                sig_cutoff))

        length = len(self.corr_rows)
        # If both sig_data and corr_data are used
        if sig_cutoff is not None and corr_cutoff is not None:
            # Go through half of matrix since keys/values repeat, keys just in reverse
            for i in range(length):
                if i % 1000 == 0:
                    logging.info('updating... on row# ' + str(i))
                for j in range(i + 1, length):
                    # Get otu's for key and then sort to cover both possibilities
                    # so won't be seperate when pushing into dict
                    sorted_otus = [self.corr_rows[i], self.corr_cols[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    # Increment frequency if
                    if self.sig_vals[i][j] <= sig_cutoff and self.corr_vals[i][j] >= corr_cutoff:
                        try:
                            del self.a_dict[key]
                        except KeyError:
                            pass

        # If no corr_cutoff is specified but sig_cutoff is; or corr_matrix doesn't exit
        elif sig_cutoff is not None:
            for i in range(length):
                for j in range(i + 1, length):
                    if i % 1000 == 0:
                        logging.info('updating... on row# ' + str(i))
                    sorted_otus = [self.sig_rows[i], self.sig_cols[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    if self.sig_vals[i][j] <= sig_cutoff:
                        try:
                            del self.a_dict[key]
                        except KeyError:
                            pass

        # if no sig_cutoff is specified but corr_cutoff is; or sig_matrix doesn't exist
        elif corr_cutoff is not None:
            for i in range(length):
                if i % 1000 == 0:
                    logging.info('updating... on row# ' + str(i))
                for j in range(i + 1, length):
                    sorted_otus = [self.corr_rows[i], self.corr_cols[j]]
                    sorted_otus.sort()
                    key = sorted_otus[0] + '<->' + sorted_otus[1]
                    if self.corr_vals[i][j] >= corr_cutoff:
                        try:
                            del self.a_dict[key]
                        except KeyError:
                            pass
        else:
            raise ValueError('ERROR: no comparing can be performed. Perhaps no corr_cutoff was specified and '
                             'significance data does not exist.')

    def _to_html(self, frequency, quantity):
        """
        Creates the html file that will display the results in a table
        """

        logging.info('to_html with frequency: {}, and quantity: {}'.format(frequency, quantity))

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
                # Values and corr and sig are averages
                if self.is_unique_search:
                    html_dict.update({key: [val[0], val[1], val[2]]})
                else:
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
        if self.sig_cutoff is not None:
            self.sig_df = pd.DataFrame(index=row_col_list, columns=row_col_list)
        self.freq_df = pd.DataFrame(index=row_col_list, columns=row_col_list)
        # Make html_str out of html_dict
        html_str = "<html>\n" \
                   "<body>\n" \
                   '<table border="2">\n' \
                   "<tr>"
        if self.is_unique_search:
            html_str += "<td>OTUs: </td><td>Correlation: </td> "
        else:
            html_str += "<td>OTUs: </td><td>Average Correlation: </td> "

        if self.sig_cutoff is not None:
            if self.is_unique_search:
                html_str += "<td>Significance: </td> "
            else:
                html_str += "<td>Average Significance: </td> "
        html_str += "<td>Frequency:</td>" \
                    "\n</tr>\n"
        counter = 0
        for key, val in html_dict.items():
            # the df part
            OTUs = key.split('<->')
            self.corr_df[OTUs[0]][OTUs[1]] = val[1]
            self.corr_df[OTUs[1]][OTUs[0]] = val[1]
            if self.sig_cutoff is not None:
                self.sig_df[OTUs[0]][OTUs[1]] = val[0]
                self.sig_df[OTUs[1]][OTUs[0]] = val[0]
            self.freq_df[OTUs[0]][OTUs[1]] = val[2]
            self.freq_df[OTUs[1]][OTUs[0]] = val[2]

            # the html part
            if counter <= 20000:
                html_str += "<tr>\n" \
                            "<td>" + key + ":</td><td>" + str(round(val[0], 5)) + "</td>"
                if self.sig_cutoff is not None:
                    html_str += "<td>" + str(round(val[1], 5)) + "</td>"
                html_str += "<td>" + str(val[2]) + " / " + str(quantity) + "</td>" \
                            + "\n</tr>\n"
            counter += 1
        html_str += "</table>\n" \
                    "</body>\n" \
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

    def run(self, MatrixIds, sig_cutoff, corr_cutoff, frequency, search_for_type, matrix_unique_to):
        pos = 1
        quantity = len(MatrixIds)
        if search_for_type == "unique":
            if matrix_unique_to is None:
                matrix_unique_to = MatrixIds[0]
            try:
                MatrixIds.remove(matrix_unique_to)
            except ValueError:
                logging.info('The matrix that was chosen for "Unique To This Matrix" was not in List of Matrices.'
                             'It will be added now.')
            MatrixIds.insert(0, matrix_unique_to)
            frequency = 1
            self._get_matrix_obj(MatrixId=matrix_unique_to, corr_cutoff=corr_cutoff, sig_cutoff=sig_cutoff)
            self._push_to_unique_dict(sig_cutoff=sig_cutoff, corr_cutoff=corr_cutoff)
            for Id in MatrixIds:
                logging.info('Analyzing matrix: {} ({} / {})'.format(Id, pos, quantity))
                pos += 1
                if Id == matrix_unique_to:
                    continue
                self._get_matrix_obj(MatrixId=Id, corr_cutoff=corr_cutoff, sig_cutoff=sig_cutoff)
                self._remove_from_unique_dict(sig_cutoff=sig_cutoff, corr_cutoff=corr_cutoff)
        if search_for_type == "union":
            frequency = 1
        if search_for_type == "intersection" or search_for_type == "union":
            if frequency > len(MatrixIds):
                raise ValueError('Frequency criteria can not have a number greater than the number of matrices '
                                 'in the list')
            for Id in MatrixIds:
                logging.info('Analyzing matrix: {} ({} / {})'.format(Id, pos, quantity))
                self._get_matrix_obj(MatrixId=Id, corr_cutoff=corr_cutoff, sig_cutoff=sig_cutoff)
                self._push_to_dict(sig_cutoff=self.sig_cutoff, corr_cutoff=self.corr_cutoff)
                pos += 1
        self._to_html(frequency=frequency, quantity=quantity)
        return {
            'html_paths': self.html_paths,
            'corr_df': self.corr_df,
            'sig_df': self.sig_df,
            'freq_df': self.freq_df
        }
