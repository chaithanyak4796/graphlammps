In_st='Dir_Lib='
Out_st=$In_st\"$(pwd)/\"
echo $Out_st
file='graphlammps/params.py'

sed -i -e 's#.*'"${In_st}"'.*#'"\\${Out_st}"'#' ${file}
find ./graphlammps -name "*.py-e" -delete

python setup.py bdist_wheel
pip install --force-reinstall dist/graphlammps-0.1.0-py3-none-any.whl
