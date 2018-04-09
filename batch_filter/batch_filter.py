#!/usr/bin/env python3

import sys

if sys.version_info < (3,5):
    print("This script requires Python >= 3.5.")
    print("Current Python version is {}.".format(sys.version.split()[0]))
    sys.exit(-1)

import time, threading, traceback, subprocess, requests
from collections import defaultdict

input_path = "all_targets.txt"
output_path = "off_targets.txt"

# For given radius c5_c10_c20 around a target we can determine if there are
# any offtargets in that radius.  There are two radii of interest.
#
offtarget_proximity = {
    "far":  "5_9_18",
    "near":  "5_9_19"
}
assert offtarget_proximity["far"] < offtarget_proximity["near"]
#
# Def: A radius c5_c10_c20 specifies maximum Hamming distances on nested
# suffixes of 5, 10, 20 bases.  A 20-mer Y is within 5_9_18 radius of 20-mer X
# if and only if the 5-char suffixes of X and Y match exactly, the 10-char
# suffixes have at most 1 positional difference, and the 20-char suffixes
# (i.e. the entire X and Y) have at most 2 positional differences.


def remove_generated_file(output_file):
    assert ' ' not in output_file
    subprocess.check_call("rm -f {}".format(output_file).split())


def fetch_with_retries(targets, c5, c10, c20, max_attempts=5, timeout=600):
    failures = 0
    while True:
        try:
            url = "http://localhost:8080/search?targets=%s&limits=%s"
            url = url % (",".join(map(str, targets)), ",".join(map(str, [c5, c10, c20])))
            return requests.get(url, timeout=timeout)
        except ConnectionResetError:
            failures += 1
            if failures > max_attempts:
                raise
            time.sleep(2 ** (failures + 5))


def run_fetch_with_retries(results, results_lock, thread_slots, arr, c5, c10, c20, i, radius):
    try:
        t = time.time()
        r = fetch_with_retries(arr, c5, c10, c20)
        with results_lock:
            results[radius].append(r)
            # print("Fetched {}-element slice at {} for {} in {:3.1f} seconds."
            #      .format(len(arr), i, radius, time.time() - t))
    except:
        traceback.print_exc()
    finally:
        thread_slots.release()


def fetch_all_offtargets(all_targets, proximity_radii):
    results = defaultdict(list)
    results_lock = threading.RLock()
    # above 8 the GO server may run out of ports or something like that
    max_threads = 8
    thread_slots = threading.Semaphore(max_threads)
    # also tried 10k, 15k, and 30k; too much variance to discern any real difference
    slice_size = 20000
    for radius in proximity_radii:
        c5, c10, c20 = tuple(int(c) for c in radius.split("_"))
        assert c5 <= 5
        assert c10 <= 10
        assert c20 <= 20
        for i in range(0, len(all_targets), slice_size):
            thread_slots.acquire()
            threading.Thread(
                target=run_fetch_with_retries,
                args=[results, results_lock, thread_slots, all_targets[i:i+slice_size], c5, c10, c20, i, radius]
            ).start()
    # wait for threads to finish
    for _ in range(max_threads):
        thread_slots.acquire()
    print("Identified all offtargets for {} x {} target spheres."
          .format(len(all_targets), len(proximity_radii)))
    return results


def output(offtargets_output, results):
    bad_targets = defaultdict(list)
    for req_text in results:
        for r in results[req_text]:
            for target_response in r.text.split('\n'):
                if target_response and target_response[0] in ('A', 'C', 'G', 'T'):
                    words = target_response.split()
                    assert len(words) == 2
                    target, boolean = words
                    assert len(target) == 20
                    assert boolean in ('false', 'true')
                    if boolean == 'true':
                        bad_targets[target].append(req_text)
                        assert len(bad_targets[target]) <= 2
    with open(offtargets_output, "w") as outf:
        for target in sorted(bad_targets.keys()):
            outf.write(target + " " + " ".join(sorted(bad_targets[target])) + "\n\n")


def read_all_targets(input_path):
    "Returns a list of the 20-mers in all_targets.txt."
    with open(input_path, "r") as f:
        targets = []
        for line in f:
            if line and line[0].isalpha():
                targets.append(line.strip())
        print("Read {} targets from {}.".format(len(targets), input_path))
        return targets


def main():
    try:
        print("Poking offtarget server.  Timeout 30 seconds.")
        fetch_with_retries(["ACGT" * 5], 5, 9, 18, max_attempts=5, timeout=30)
        print("Offtarget server is alive.")
    except:
        traceback.print_exc()
        print("*********************************************************************************")
        print("***   Did you forget to start the offtarget server?  Did it finish loading?   ***")
        print("***   Please follow the instructions in README.TXT and try again.             ***")
        print("*********************************************************************************")
        sys.exit(-1)
    t = time.time()
    remove_generated_file(output_path)
    all_targets = read_all_targets(input_path)
    results = fetch_all_offtargets(all_targets, offtarget_proximity.values())
    output(output_path, results)
    print("Completed filter_offtarget in {:3.1f} seconds.".format(time.time() - t))


if __name__ == "__main__":
    retcode = main()
    sys.exit(retcode)
