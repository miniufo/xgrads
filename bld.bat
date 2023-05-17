echo "Building xgrads to conda (win)"
"%PYTHON%" setup.py install
if errorlevel 1 exit 1
