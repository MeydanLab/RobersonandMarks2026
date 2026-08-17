import os
import sys

import tools_m


if __name__ == "__main__":
	if len(sys.argv) > 1:
		infile=sys.argv[1]
	else:
		here=os.path.dirname(os.path.abspath(__file__))
		infile=os.path.join(here,"inputfiles","metagene_transcript.txt")
	tools_m.metagene_transcript_m_wrapper(infile)
