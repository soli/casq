"""Various utilities."""
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
    process = subprocess.run(command.split(), capture_output=True, check=True)
    result = json.loads(process.stdout)["validation-results"]
    if "no-errors" in result:
        return "OK"
    else:
        return json.dumps(result["problem"], indent=2)
