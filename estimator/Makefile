SAGE ?= sage
OUT ?= LoTRS-Estimate-Output-N100T50.generated.txt

.PHONY: all estimate reference clean

all: estimate

estimate:
	$(SAGE) -c "exec(open('lotrs_estimate.py').read()); main()"

reference:
	$(SAGE) -c "exec(open('lotrs_estimate.py').read()); main()" | tee $(OUT)

clean:
	$(RM) -f *.pyc *.cprof */*.pyc *.rsp *.log
	$(RM) -rf __pycache__ */__pycache__

