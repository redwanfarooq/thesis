#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>

using namespace std;
using namespace Rcpp;

// Progress display function with Jupyter notebook compatibility
inline void show_progress(int current, int total, bool& completed) {
    if (completed) return;
    
    // Use cat() and flush.console() for immediate display in Jupyter notebooks
    Function cat("cat");
    Function flush_console("flush.console");
    
    if (total > 0) {
        // Known total - show percentage progress bar
        float progress = (float)current / total;
        int bar_width = 40;
        int pos = (int)(bar_width * progress);
        int percent = (int)(progress * 100);
        
        // Build progress bar string
        string bar = "\r[";
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) bar += "=";
            else if (i == pos) bar += ">";
            else bar += " ";
        }
        bar += "] " + to_string(percent) + "% (" + 
               to_string(current) + "/" + to_string(total) + ")";
        
        cat(bar);
        flush_console();
        
        // Mark completed when 100% is reached
        if (current >= total) {
            cat("\n");
            flush_console();
            completed = true;
        }
    } else {
        // Unknown total - show simple counter
        string progress_msg = "\rProcessed " + to_string(current) + " rows...";
        cat(progress_msg);
        flush_console();
    }
}

// High-performance line parsing and statistics calculation
inline pair<double, double> parse_line_and_calculate_stats(const char* buffer, size_t length, char delimiter = '\t') {
    const char* ptr = buffer;
    const char* end = buffer + length;
    
    // Use register variables for frequently accessed data
    register const char* p = ptr;
    register double s = 0.0, ss = 0.0;
    register int c = 0;
    
    while (p < end) {
        // Skip delimiters
        while (p < end && *p == delimiter) ++p;
        if (p >= end) break;
        
        // Fast number parsing
        char* next_ptr;
        double val = strtod(p, &next_ptr);
        
        // Update statistics if parsing succeeded
        if (next_ptr > p) {
            s += val;
            ss += val * val;
            ++c;
        }
        
        p = next_ptr;
    }
    
    if (c == 0) {
        return make_pair(R_NaN, R_NaN);
    }
    
    // Calculate mean and standard deviation
    double inv_n = 1.0 / c;
    double mean = s * inv_n;
    
    if (c == 1) {
        return make_pair(mean, 0.0);
    }
    
    // Numerically stable variance calculation
    double variance = (ss - s * mean) / (c - 1);
    double std_dev = variance > 0.0 ? sqrt(variance) : 0.0;
    
    return make_pair(mean, std_dev);
}

//' Calculate row-wise statistics for TSV files
//' 
//' Efficiently processes TSV files (including GZIP compressed) line by line
//' to calculate mean and standard deviation for each row.
//' 
//' @param file_path Path to the TSV file (.tsv or .tsv.gz)
//' @param max_rows Maximum number of rows to process (0 = all rows)
//' @return DataFrame with columns: row_index, mean, sd
//' @export
// [[Rcpp::export]]
DataFrame tsv_row_stats(string file_path, int max_rows = 0) {
    
    // Initialize result vectors
    vector<double> means;
    vector<double> std_devs;
    
    // Pre-allocate vectors for better performance
    if (max_rows > 0) {
        means.reserve(max_rows);
        std_devs.reserve(max_rows);
    } else {
        // Conservative estimate for unknown size
        means.reserve(1000);
        std_devs.reserve(1000);
    }
    
    string line;
    line.reserve(2048); // Pre-allocate for typical line lengths
    int processed_rows = 0;
    bool progress_completed = false;
    
    // Check if file is GZIP compressed
    bool is_gzip = (file_path.substr(file_path.find_last_of(".") + 1) == "gz");
    
    // Setup progress tracking
    int total_lines = 0;
    int expected_rows = 0;
    
    if (max_rows == 0) {
        Function message("message");
        message("Processing all rows (progress will be shown every 1000 rows)...");
        total_lines = -1; // Unknown total
        expected_rows = -1;
    } else {
        expected_rows = max_rows;
        total_lines = max_rows;
    }
    
    // Process file based on compression type
    if (is_gzip) {
        // Handle GZIP compressed files using gunzip pipe
        string command = "gunzip -c '" + file_path + "'";
        FILE* pipe = popen(command.c_str(), "r");
        if (!pipe) {
            stop("Failed to open compressed file: " + file_path);
        }
        
        // Variables for efficient line reading from pipe
        size_t buffer_size = 0;
        char* buffer = nullptr;
        ssize_t line_length;
        
        try {
            while ((line_length = getline(&buffer, &buffer_size, pipe)) != -1) {
                // Remove trailing newline
                if (line_length > 0 && buffer[line_length-1] == '\n') {
                    buffer[line_length-1] = '\0';
                    line_length--;
                }
                
                // Skip empty lines
                if (line_length == 0) continue;
                
                // Check row limit
                if (max_rows > 0 && processed_rows >= max_rows) break;
                
                // Calculate statistics for this row
                pair<double, double> stats = parse_line_and_calculate_stats(buffer, line_length);
                
                if (!isnan(stats.first)) {
                    means.push_back(stats.first);
                    std_devs.push_back(stats.second);
                    processed_rows++;
                    
                    // Progress updates
                    int progress_interval = (total_lines > 0) ? 50 : 1000;
                    if (processed_rows % progress_interval == 0 || (total_lines > 0 && processed_rows >= min(expected_rows, total_lines))) {
                        Rcpp::checkUserInterrupt();
                        show_progress(processed_rows, (total_lines > 0) ? min(expected_rows, total_lines) : -1, progress_completed);
                    }
                }
            }
            
            if (buffer) free(buffer);
            
        } catch (...) {
            if (buffer) free(buffer);
            pclose(pipe);
            throw;
        }
        
        pclose(pipe);
        
    } else {
        // Handle uncompressed files
        ifstream infile(file_path);
        if (!infile.is_open()) {
            stop("Failed to open file: " + file_path);
        }
        
        while (getline(infile, line)) {
            // Skip empty lines
            if (line.empty()) continue;
            
            // Check row limit
            if (max_rows > 0 && processed_rows >= max_rows) break;
            
            // Calculate statistics for this row
            pair<double, double> stats = parse_line_and_calculate_stats(line.c_str(), line.length());
            
            if (!isnan(stats.first)) {
                means.push_back(stats.first);
                std_devs.push_back(stats.second);
                processed_rows++;
                
                // Progress updates
                int progress_interval = (total_lines > 0) ? 50 : 1000;
                if (processed_rows % progress_interval == 0 || (total_lines > 0 && processed_rows >= min(expected_rows, total_lines))) {
                    Rcpp::checkUserInterrupt();
                    show_progress(processed_rows, (total_lines > 0) ? min(expected_rows, total_lines) : -1, progress_completed);
                }
            }
        }
        
        infile.close();
    }
    
    // Show completion message for unknown totals
    if (total_lines <= 0 && processed_rows > 0) {
        Function cat("cat");
        Function flush_console("flush.console");
        cat("\rCompleted processing " + to_string(processed_rows) + " rows.\n");
        flush_console();
    }
    
    // Create and return DataFrame
    return DataFrame::create(
        Named("mean") = means,
        Named("sd") = std_devs,
        Named("stringsAsFactors") = false
    );
}
