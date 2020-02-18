/*
A KBase module: SignificantInteractions
*/

module SignificantInteractions {
    typedef structure {
        string MatrixIds;
        float cutoff;
        int frequency;
        string corr_matrix_name;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_SignificantInteractions(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
