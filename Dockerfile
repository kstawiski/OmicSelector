FROM kstawiski/omicselector-gpu
# speed up for ainize

EXPOSE 8888
EXPOSE 80
EXPOSE 443

ENTRYPOINT ["/entrypoint.sh"]