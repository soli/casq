"""Various utilities.

Copyright (C) 2019 Sylvain Soliman

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import json
import subprocess


def validate(filename: str) -> str:
    """Validate an SBML file using the online validator API.

    Unit consistency verification is off
    """
    url = "http://sbml.org/validator/"
    command = "curl -s -F file=@{filename} -F output=json -F offcheck=u {url}".format(
        filename=filename, url=url
    )
    process = subprocess.run(command.split(), stdout=subprocess.PIPE, check=True)
    result = json.loads(process.stdout.decode("utf-8"))["validation-results"]
    if "no-errors" in result:
        return "OK"
    # else:
    return json.dumps(result["problem"], indent=2)
