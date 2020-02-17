/*
A KBase module: SignificantInteractions
*/

module SignificantInteractions {
    typedef structure {
        string MatrixIds;
        float cutoff;
        int frequency;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_SignificantInteractions(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
