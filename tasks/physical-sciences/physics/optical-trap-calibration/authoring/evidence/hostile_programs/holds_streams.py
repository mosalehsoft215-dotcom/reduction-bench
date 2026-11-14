import os, time
if not os.fork():
    os.setsid()
    time.sleep(3000)
    os._exit(0)
