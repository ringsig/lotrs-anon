clean:
	cd lotrs-rs && cargo clean
	cd lotrs-py && $(MAKE) clean
	cd estimator && $(MAKE) clean

