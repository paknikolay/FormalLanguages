#include <iostream>
#include <vector>
#include <string>
using std::string;
using std::vector;
using std::cout;
using std::cin;

class KeeperOfDeduciblePrefixes {
private:
    string *word_;
    vector<vector<int>> deducibleSubsequencesFrefixes; //deducible prefixes
    vector<vector<int>> deducibleSubsequences; //deducible subsequences

public:

    explicit KeeperOfDeduciblePrefixes(string *word) : word_(word) {
        for (int i = 0; i <= word_->length(); ++i) {//initialization
            deducibleSubsequencesFrefixes.emplace_back(vector<int>());
            deducibleSubsequences.emplace_back(vector<int>());
            deducibleSubsequencesFrefixes[i].assign(word_->length() + 1, false);
            deducibleSubsequences[i].assign(word_->length() + 1, false);
        }
    }

    KeeperOfDeduciblePrefixes(string *word, char symbol) : word_(word) {
        for (int i = 0; i <= word_->length(); ++i) {//initialization
            deducibleSubsequencesFrefixes.emplace_back(vector<int>());
            deducibleSubsequencesFrefixes[i].assign(word_->length() + 1, false);
            deducibleSubsequences.emplace_back(vector<int>());
            deducibleSubsequences[i].assign(word_->length() + 1, false);
        }

        if (symbol == '1') {
            for (int j = 0; j < word_->length() + 1; ++j) {
                deducibleSubsequencesFrefixes[j][0] = true;
                deducibleSubsequences[j][0] = true;
            }
        } else {
            for (int i = 0; i < word_->length(); ++i) {
                if ((*word_)[i] == symbol) {
                    deducibleSubsequencesFrefixes[i][1] = true;
                    deducibleSubsequences[i][1] = true;
                }
            }
        }
    }

    KeeperOfDeduciblePrefixes(const KeeperOfDeduciblePrefixes &other) : word_(other.word_) {
        for (int i = 0; i <= word_->length(); ++i) {
            deducibleSubsequencesFrefixes.emplace_back(vector<int>());
            deducibleSubsequences.emplace_back(vector<int>());
            for (int j = 0; j <= word_->length(); ++j) {
                deducibleSubsequencesFrefixes[i].emplace_back(deducibleSubsequencesFrefixes[i][j]);
                deducibleSubsequences[i].emplace_back(deducibleSubsequences[i][j]);
            }
        }
    }

    KeeperOfDeduciblePrefixes &operator=(const KeeperOfDeduciblePrefixes &other) {
        word_ = other.word_;
        for (int i = 0; i <= word_->length(); ++i) {
            deducibleSubsequencesFrefixes.emplace_back(vector<int>());
            deducibleSubsequences.emplace_back(vector<int>());
            for (int j = 0; j <= word_->length(); ++j) {
                deducibleSubsequences[i].emplace_back(deducibleSubsequences[i][j]);
                deducibleSubsequencesFrefixes[i].emplace_back(deducibleSubsequencesFrefixes[i][j]);
            }
        }
    }

    KeeperOfDeduciblePrefixes operator+(const KeeperOfDeduciblePrefixes &other) const { //union deducible prefixes
        KeeperOfDeduciblePrefixes tmp(word_);
        for (int i = 0; i <= word_->length(); ++i) {
            for (int j = 0; j <= word_->length(); ++j) {
                tmp.deducibleSubsequencesFrefixes[i][j] =
                        deducibleSubsequencesFrefixes[i][j] | other.deducibleSubsequencesFrefixes[i][j];
                tmp.deducibleSubsequences[i][j] = deducibleSubsequences[i][j] | other.deducibleSubsequences[i][j];
            }
        }
        return tmp;
    }

    //* means concatenation when  a*b
    KeeperOfDeduciblePrefixes operator*(const KeeperOfDeduciblePrefixes &other) const {
        KeeperOfDeduciblePrefixes tmp(word_);
        for (int i = 0; i <= word_->length(); ++i) {
            for (int j = 0; j <= word_->length(); ++j) {
                tmp.deducibleSubsequencesFrefixes[i][j] = false;
                tmp.deducibleSubsequences[i][j] = false;
            }
        }

        for (int i = 0; i < word_->length(); ++i) {
            for (int j = 0; j <= word_->length(); ++j) {
                if (deducibleSubsequences[i][j]) //if there is such subsequence
                    for (int k = 0; k <= word_->length(); ++k) { //trying to prolong current subsequence
                        if (other.deducibleSubsequences[i + j][k]) {
                            tmp.deducibleSubsequences[i][j + k] = true;
                            tmp.deducibleSubsequencesFrefixes[i][j + k] = true;
                        }
                    }
            }
        }
        //adding prefix in the end
        for (int i = 0; i <= word_->length(); ++i) {
            for (int j = 0; j <= word_->length(); ++j) {
                if (deducibleSubsequences[i][j]) { //if there is such subsequence
                    for (int k = 0; k <= word_->length(); ++k) { //trying to prolong current subsequence
                        if (other.deducibleSubsequencesFrefixes[i + j][k])
                            tmp.deducibleSubsequencesFrefixes[i][j + k] = true;
                    }
                }
                if (deducibleSubsequencesFrefixes[i][j])
                    tmp.deducibleSubsequencesFrefixes[i][j] = true;
            }
        }

        return tmp;
    }

    //* means Kleene star when *a
    KeeperOfDeduciblePrefixes operator*() const {
        KeeperOfDeduciblePrefixes tmp(word_);
        for (int i = 0; i <= word_->length(); ++i) {
            for (int j = 0; j <= word_->length(); ++j) {
                tmp.deducibleSubsequencesFrefixes[i][j] = deducibleSubsequences[i][j];
                tmp.deducibleSubsequences[i][j] = deducibleSubsequences[i][j];
            }
            tmp.deducibleSubsequencesFrefixes[i][0] = true;
            tmp.deducibleSubsequences[i][0] = true;
        }

        for (int i = 0; i <= word_->length(); ++i) {
            for (int j = 0; j <= word_->length(); ++j) {
                if (tmp.deducibleSubsequences[i][j]) //if there is such subsequence
                    for (int k = 0; k <= word_->length(); ++k) { //trying to prolong current subsequence
                        if (tmp.deducibleSubsequences[i + j][k]) {
                            tmp.deducibleSubsequencesFrefixes[i][j + k] = true;
                            tmp.deducibleSubsequences[i][j + k] = true;

                        }
                    }
            }
        }
        //adding prefix in the end
        for (int i = 0; i <= word_->length(); ++i) {
            for (int j = 0; j <= word_->length(); ++j) {
                if (tmp.deducibleSubsequences[i][j]) //if there is such subsequence
                    for (int k = 0; k <= word_->length(); ++k) { //trying to prolong current subsequence
                        if (deducibleSubsequencesFrefixes[i + j][k]) tmp.deducibleSubsequencesFrefixes[i][j + k] = true;
                    }
            }
        }

        return tmp;
    }

    int getMaxPrefix() {
        int maxLength = 0;
        for (int i = 0; i <= word_->length(); ++i) {
            if (deducibleSubsequencesFrefixes[0][i]) maxLength = i;
        }
        return maxLength;
    }
};


class MaximumPrefixFinder {
private:
    string *word_, *regularExpression_;
    vector<KeeperOfDeduciblePrefixes *> links_;
    vector<char> stack_;

    bool checkIfIsValid() {//check valid symbols
        for (int i = 0; i < word_->length(); ++i) {
            char symbol = (*word_)[i];
            if (!((symbol >= 'a' && symbol <= 'c') || symbol == '1'))return false;
        }

        for (int i = 0; i < regularExpression_->length(); ++i) {
            char symbol = (*regularExpression_)[i];
            if (!((symbol >= 'a' && symbol <= 'c') || symbol == '1' || symbol == '+' || symbol == '*' ||
                  symbol == '.'))
                return false;
        }
        return true;
    }

public:

    explicit MaximumPrefixFinder(string *regularExpression, string *word) : word_(word),
                                                                            regularExpression_(regularExpression) {
    }

    int findMaxLengthOfPrefix() {
        checkIfIsValid();
        int iter = 0; //pointer to position in the regular expression
        while (iter < regularExpression_->length()) {
            char currentSymbol = (*regularExpression_)[iter];
            if (currentSymbol >= 'a' && currentSymbol <= 'z' || currentSymbol == '1') {
                stack_.push_back(currentSymbol);
                links_.push_back(nullptr);
            } else {
                if (stack_.empty()) return -1;//error
                if (currentSymbol == '*') {
                    if (links_.back() == nullptr) {
                        links_.back() = new KeeperOfDeduciblePrefixes(word_, stack_.back());
                    }
                    KeeperOfDeduciblePrefixes *tmp = new KeeperOfDeduciblePrefixes(*(*links_.back()));//applying *

                    //deleting unnecessary links and symbols in stack_
                    delete links_.back();
                    stack_.back() = 'r'; // regular expression
                    links_.back() = tmp;
                } else {
                    if (stack_.size() < 2) return -1; //not enough arguments

                    if (links_.back() == nullptr) {
                        links_.back() = new KeeperOfDeduciblePrefixes(word_, stack_.back());
                    }

                    if (links_[links_.size() - 2] == nullptr) {
                        links_[links_.size() - 2] = new KeeperOfDeduciblePrefixes(word_, stack_[stack_.size() - 2]);
                    }

                    KeeperOfDeduciblePrefixes *tmp;

                    if (currentSymbol == '+') {
                        tmp = new KeeperOfDeduciblePrefixes(*links_[links_.size() - 2] + *links_.back());//applying +
                    } else {
                        tmp = new KeeperOfDeduciblePrefixes(*links_[links_.size() - 2] * *links_.back());//applying *
                    }

                    //deleting unnecessary links and symbols in stack_
                    delete links_.back();
                    links_.pop_back();
                    delete links_.back();
                    links_.back() = tmp;

                    stack_.pop_back();
                    stack_.back() = 'r'; //regular expression
                }
            }
            ++iter;
        }
        if (stack_.size() != 1)
            return -1;// error

        int maxLength = 0;
        return links_.back()->getMaxPrefix();
    }
};

int main() {

    string regexpr, word;
    cin >> regexpr >> word;
    MaximumPrefixFinder maximumPrefixFinder(&regexpr, &word);
    cout << maximumPrefixFinder.findMaxLengthOfPrefix();
    return 0;
}