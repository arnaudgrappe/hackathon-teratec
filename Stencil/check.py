import sys

def process_file(filename):
	with open(filename) as f:
		lines = f.readlines()

	for i in range(len(lines)):
		lines[i] = lines[i].split()[1:]
		for j in range(len(lines[i])):
			lines[i][j] = float(lines[i][j])
	
	return lines


if __name__ == "__main__":

	# check arguments
	if (len(sys.argv) <= 2):
		print("Usage: python3 check.py <job.out> <ref.out>")
		exit()

	job_file = sys.argv[1]
	ref_file = sys.argv[2]
	ref = process_file(ref_file)
	job = process_file(job_file)

	# check if dimensions match
	if (ref[0][7:10] != job[0][7:10]):
		print("ERROR: mismatched dimensions.")
		print(f"{job_file}: {int(job[0][7])} {int(job[0][8])} {int(job[0][9])}")
		print(f"{ref_file}: {int(ref[0][7])} {int(ref[0][8])} {int(ref[0][9])}")
		exit()

	# check if both files have as many iterations
	iterations = min(len(job), len(ref))
	if (len(ref) != len(job)):
		print("WARNING: not the same number of iterations.")
		print(f"  {job_file}: {len(job)}")
		print(f"  {ref_file}: {len(ref)}")
		print(f"  The first {iterations} iterations will be compared.\n")

	dimx, dimy, dimz = [int(dim) for dim in ref[0][7:10]]
	print(f"Config: {dimx} {dimy} {dimz} ({iterations} iterations)")

	mean_error    = 0
	max_error     = 0
	mean_time     = 0
	mean_time_ref = 0
	mean_speedup  = 0
	for i in range(iterations):
		for j in range(5):
			error = abs(job[i][j] - ref[i][j])
			if error > max_error:
				max_error = error
			mean_error += error
		mean_time += job[i][5]
		mean_time_ref += ref[i][5]

	mean_time /= iterations
	mean_time_ref /= iterations
	mean_error /= iterations * 5
	mean_speedup = mean_time_ref / mean_time

	print("--------------------")
	print(f"Mean error  {mean_error:.2e}")
	print(f"Max error   {max_error:.2e}")
	print("--------------------")
	print(f"Mean time   {mean_time / 1000:.3f} ms")
	print(f"Speedup     {mean_speedup:.1f}")