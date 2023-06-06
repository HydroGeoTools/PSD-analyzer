FROM python:3.10

# Set up working dir and copy app
RUN mkdir /app
WORKDIR /app
ADD . /app

# Set up Python app environment
RUN pip install --no-cache-dir -r requirements.txt

# Launch
CMD ["gunicorn", "wsgi:application"]
