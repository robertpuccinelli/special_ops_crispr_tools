#include "catch.hpp"

#include <unistd.h>

#include <cstdlib>
#include <vector>
#include <iostream>
#include <assert.h>
#include <errno.h>
#include <fcntl.h>

#include "../crispr_sites.hpp"

using namespace std;

// unit tests for scan_stdin()

// forward declarations we need
void scan_stdin(bool output_counts);

void decode(char* buf, const int len, const int64_t code);
template<int len> int64_t encode(const char* buf);
int64_t complement(const int64_t code);
void init_encoding();

char random_base(void) {
    switch (rand() % 4) {
    case 0:
	return 'A';
    case 1:
	return 'C';
    case 2:
	return 'G';
    default:
	return 'T';
    }
}

void random_kmer(char* output) {
    for (int i = 0; i < k; i++) {
	output[i] = random_base();
    }
}

void random_sequence_no_pam(char* output, int len) {
    for (int i = 0; i < len; i++) {
	if ((i > 0) && (output[i - 1] == 'C' || output[i - 1] == 'G')) {
	    output[i] = 'A';
	} else {
	    output[i] = random_base();
	}
    }
}

TEST_CASE( "scan_stdin correctly finds crispr sites", "[scan_stdin]" ) {
    init_encoding();
    
    constexpr auto READ_PIPE = 0;
    constexpr auto WRITE_PIPE = 1;

    // save original stdin and stdout
    int original_stdin;
    int original_stdout;
    int original_stderr;

    int err;

    original_stdin = dup(STDIN_FILENO);

    assert(original_stdin != -1);
    
    original_stdout = dup(STDOUT_FILENO);

    assert(original_stdout != -1);

    original_stderr = dup(STDERR_FILENO);

    assert(original_stderr != -1);
    
    // hook up stdin and stdout to pipes we control
    int stdin_pipe[2];
    int stdout_pipe[2];
    int stderr_pipe[2];
    
    err = pipe(stdin_pipe);

    assert(err != -1);

    err = pipe(stdout_pipe);

    assert(err != -1);

    err = pipe(stderr_pipe);

    assert(err != -1);

    err = dup2(stdin_pipe[READ_PIPE], STDIN_FILENO);

    assert(err != -1);

    err = dup2(stdout_pipe[WRITE_PIPE], STDOUT_FILENO);

    assert(err != -1);

    err = dup2(stderr_pipe[WRITE_PIPE], STDERR_FILENO);

    assert(err != -1);
    
    vector<string> expected_crispr_sites;

    constexpr auto INPUT_SIZE = 1024;
    char input[INPUT_SIZE];

    // We need enough room in the INPUT to place test cases at various positions
    assert(INPUT_SIZE > 10*k);
    
    // Initialize input to a random sequence containing no crispr sites
    random_sequence_no_pam(input, INPUT_SIZE);

    SECTION("detect nothing if given no crispr sites") {
    	// noop, just make sure 
    }
    SECTION("expected crispr sites detected") {
    	for (int i = k - 1; i < INPUT_SIZE - k; i++) {
    	    if (rand() % 10 == 0) {
    		input[i + 1] = 'T';
    		input[i] = 'G';
    		input[i - 1] = 'G';
    		input[i - 2] = 'A';

    		char crispr_site[k - 2];
    		crispr_site[k - 3] = 0;

    		for (int j = 0; j < k - 3; j++) {
    		    crispr_site[j] = input[i - (k - 1) + j];
    		}
    		expected_crispr_sites.push_back(string(crispr_site));

    		// spread out crispr sites
    		i += k + 1;
    	    }
    	}
    }
    SECTION("crispr sites not detected across separators") {
    	int separator_position = rand() % (INPUT_SIZE - 5*k) + 2*k;

    	input[separator_position] = '>';
    	input[separator_position + k] = '\n';

    	input[separator_position + k + 3] = 'G';
    	input[separator_position + k + 2] = 'G';	

    	input[separator_position - 2] = 'C';
    	input[separator_position - 1] = 'C';	
    }
    SECTION("expected crispr sites detected with separators") {
	for (int i = k; i < INPUT_SIZE - 3*k - 2; i++) {
	    if (rand() % 10 == 0) {
		input[i] = 'T';
		input[i - 1] = 'G';
		input[i - 2] = 'G';
		input[i - 3] = 'A';

		char crispr_site[k - 2];
		crispr_site[k - 3] = 0;

		for (int j = 0; j < k - 3; j++) {
		    crispr_site[j] = input[i - k + j];
		}
		cerr << i << " " << crispr_site << endl;
		expected_crispr_sites.push_back(string(crispr_site));

		// spread out crispr sites
		input[i + 2] = '>';
		input[i + 2*k] = '\n';
		i += 3*k + 2;
	    }
	}
    }
    
    // Run the program!
    write(stdin_pipe[WRITE_PIPE], input, strlen(input));
    close(stdin_pipe[WRITE_PIPE]);

    scan_stdin(false);

    fflush(stdout);

    char buffer[20];

    fcntl(stdout_pipe[READ_PIPE], F_SETFL, O_NONBLOCK);

    vector<string> detected_crispr_sites;
    
    // Read redirected stdin until the stream is empty        
    while (true) {
	char crispr_site[k - 2];
	crispr_site[k - 3] = 0;
	
        const ssize_t bytes_read = read(stdout_pipe[READ_PIPE], crispr_site, k - 3);

	if (bytes_read == -1) {
	    if (errno == EAGAIN || errno == EWOULDBLOCK) {
		// we finished reading everything from redirected stdin
		break;
	    }

	    // unexpected error encoutered
	    assert(false);
	}

        // end of file	
        if (bytes_read == 0) {
            break;
        }

	detected_crispr_sites.push_back(string(crispr_site));

	// Read & discard newline
	read(stdout_pipe[READ_PIPE], crispr_site, 1);
    }
    
    close(stdout_pipe[READ_PIPE]);


    // Restore original stdin/stdout
    dup2(original_stdin, STDIN_FILENO);
    dup2(original_stdout, STDOUT_FILENO);
    dup2(original_stderr, STDERR_FILENO);
    
    REQUIRE(detected_crispr_sites.size() == expected_crispr_sites.size());
    
    for (int i = 0; i < detected_crispr_sites.size(); i++) {
	bool valid_detection = false;
	for (int j = 0; j < expected_crispr_sites.size(); j++) {
	    if (expected_crispr_sites[j] == detected_crispr_sites[i]) {
		valid_detection = true;
		break;
	    }
	}
	REQUIRE(valid_detection == true);
    }
}
