#include "headers/compat.h"

// symbols are automatically created and stored in the table
// this way, when a (for instance) letter shows up more than once,
// GiNaC will still see interpret it as the same symbol object
g::symtab table;
parser reader(table);

g::matrix toMatrix(string str) {
    compress_whitespace_and_newlines(str);
    string meat = str.substr(1, str.size() - 2);
    //std::cout << meat << std::endl;
    vector<string> rows = split(meat, ";");

    // Make it into a 2d string vector ....
    // Can probably combine this with the next step.
    vector<vector<string>> mat = {};
    for (int i = 0; i < rows.size(); i++) {
        string rowString = rows[i];
        trim(rowString);
        // If there are commas, use commas as the delimiter
        string delimiter = (rowString.find(',') != std::string::npos) ? "," : " ";
        vector<string> row = split(rowString, delimiter);
        mat.push_back(row);
    }

    int r = mat.size();
    int c = mat[0].size();
    list<ex> els = {};
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            els.push_back(reader(mat[i][j]));
        }
    }
    lst myLst(els);

    g::matrix result = g::matrix(r, c, myLst);

    return result;
}

string toString(g::matrix a) {
    string result = "";
    for (int i = 0; i < a.rows(); i++) {
        if (i != 0) result += ";";
        for (int j = 0; j < a.cols(); j++) {
            if (j != 0) result += " ";
            std::ostringstream s;
            s << g::dflt << a(i, j);
            result += s.str();
        }
    }
    return "[" + result + "]";
}

vector<g::matrix> toMatrixSequence(string str) {
    compress_whitespace_and_newlines(str);
    vector<string> matrixStrings = split(str, "\n");

    vector<g::matrix> seq = vector<g::matrix>();
    for (int i = 0; i < matrixStrings.size(); i++) {
        trim(matrixStrings[i]);
        g::matrix matrix = toMatrix(matrixStrings[i]);
        seq.push_back(matrix);
    }

    return seq;
}

string toString(vector<g::matrix> seq) {
    string seqString = "";
    for (int i = 0; i < seq.size(); i++) {
        if (i != 0) {
            seqString += "\n";
        }
        seqString += toString(seq[i]);
    }

    return seqString;
}

lie_algebra* toLieAlgebra(string str) {
    compress_whitespace_and_newlines(str);
    vector<g::matrix> generators = toMatrixSequence(str);
    lie_algebra* l = new lie_algebra(generators, false);
    return l;
}

string toString(lie_algebra* l) {
    vector<g::matrix> basis = l->get_basis();
    return toString(basis);
}

vector<lie_algebra*> toLieAlgebraSequence(string str) {
    compress_whitespace_and_newlines(str);
    vector<string> algebraStrings = split(str, "\n@\n");

    vector<lie_algebra*> seq = vector<lie_algebra*>();
    for (int i = 0; i < algebraStrings.size(); i++) {
        trim(algebraStrings[i]);
        lie_algebra* algebra = toLieAlgebra(algebraStrings[i]);
        seq.push_back(algebra);
    }

    return seq;
}

string toString(vector<lie_algebra*> seq) {
    string seqString = "";
    for (int i = 0; i < seq.size(); i++) {
        if (i != 0) {
            seqString += "\n@\n";
        }
        seqString += toString(seq[i]);
    }

    return seqString;
}

g::matrix toMatrix(const char* s) {
    string str(s);
    return toMatrix(str);
}

const char* toCharArray(g::matrix a) {
    string str = toString(a);
    return toCharArray(str);
}

vector<g::matrix> toMatrixSequence(const char* s) {
    string str(s);
    return toMatrixSequence(str);
}

const char* toCharArray(vector<g::matrix> b) {
    string str = toString(b);
    return toCharArray(str);
}

lie_algebra* toLieAlgebra(const char* s) {
    string str(s);
    return toLieAlgebra(str);
}

const char* toCharArray(lie_algebra* l) {
    string str = toString(l);
    return toCharArray(str);
}

vector<lie_algebra*> toLieAlgebraSequence(const char* s) {
    string str(s);
    return toLieAlgebraSequence(str);
}

const char* toCharArray(vector<lie_algebra*> b) {
    string str = toString(b);
    return toCharArray(str);
}

// UTIL
const char* toCharArray(string str) {
    // declaring character array (+1 for null terminator)
    char* result = new char[str.length() + 1];
    // copying the contents of the
    // string to char array
    strcpy(result, str.c_str());
    return result;
}

string join(vector<string> strings, string joiner) {
    string result = "";
    for (int i = 0; i < strings.size(); i++) {
        if (i != 0) {
            result += joiner;
        }
        result += strings[i];
    }
    return result;
}

vector<string> split(string target, string delimiter) {
    vector<string> components;
    if (!target.empty()) {
        size_t start = 0;
        do {
            const size_t index = target.find(delimiter, start);
            if (index == string::npos) {
                break;
            }
            const size_t length = index - start;
            components.push_back(target.substr(start, length));
            start += (length + delimiter.size());
        } while (true);
        components.push_back(target.substr(start));
    }

    return components;
}

string replace(string target, string oldseg, string newseg) {
    vector<string> components = split(target, oldseg);
    return join(components, newseg);
}

void ltrim(string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

void rtrim(string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void trim(string &s) {
    rtrim(s);
    ltrim(s);
}

bool BothAreWhitespace(char lhs, char rhs) { return std::isspace(lhs) && std::isspace(rhs); }
bool BothAreNewlines(char lhs, char rhs) { return (lhs == rhs) && (lhs == '\n'); }

void compress_whitespace_and_newlines(string &s) {
    std::string::iterator new_end = std::unique(s.begin(), s.end(), BothAreWhitespace);
    s.erase(new_end, s.end());
    new_end = std::unique(s.begin(), s.end(), BothAreNewlines);
    s.erase(new_end, s.end());
}
